# useful functions for the lipid kinetic prototype
import os
import numpy as np
import pylab as plt
from scipy.optimize import minimize
from molmass import Formula
from mass_spec_utils.adduct_calculator.adduct_rules import AdductTransformer, ParsingAdductTransformer
import xlsxwriter


def min_func(x, *args):
    # x = k: the rate constant
    # *args
    p = args
    if p[0] == True:
        # fix ends is true - fix the end points from the data - not recommended
        t = p[1]
        data = p[2]
        a0 = p[3]
        ai = p[4]

        k = x[0]
    else:
        t = p[1]
        data = p[2]

        k = x[0]
        a0 = x[1]
        ai = x[2]

    preds = a0 + (ai-a0)*(1-np.exp(-k*t))
    return ((preds-data)**2).sum()


def fit(times, data_mat, fix_ends=True, make_plot=True, options={}, method='L-BFGS-B'):

    t = np.array(times)

    a0 = data_mat[0, 0]
    ai = data_mat[-1, 0]

    k = 0.05

    if fix_ends:
        args = (fix_ends, t, data_mat[:, 0], a0, ai)
        zinit = [k]
    else:
        args = (fix_ends, t, data_mat[:, 0])
        zinit = [k, a0, ai]

    end = minimize(min_func, zinit, args=args, method=method)

    if fix_ends:
        k = end['x'][0]
    else:
        k = end['x'][0]
        a0 = end['x'][1]
        ai = end['x'][2]

    if make_plot:
        plt.figure()
        plt.plot(t, data_mat[:, 0], 'ro')
        plt.plot(t, a0 + (ai-a0)*(1-np.exp(-k*t)))

    return (k, a0, ai)


def get_iso_intense(mzml_file_obj, target_rt_range, formula, adduct_type, mz_tol=(5, 'ppm'), scan_delta=2, max_iso_n=5):

    relevant_scans = list(filter(lambda x: x.rt_in_seconds >= target_rt_range[0] and
                                 x.rt_in_seconds <= target_rt_range[1], mzml_file_obj.scans))
    spectrum = Formula(formula).spectrum()
    if adduct_type == '[M-H+FA]-':
        adduct_type = '[M-H+CH2O2]-'
    target_mz = [AdductTransformer().mass2ion(x[0], adduct_type)
                 for x in spectrum.values()]
    target_mz.sort()
    isos = []
    t = target_mz[0]
    max_i = 0
    max_mz = -1
    max_idx = None
    max_rt = None
    max_scan_no = None

    if mz_tol[1] == 'ppm':
        mz_tol_abs = t*mz_tol[0]/1e6
    else:
        mz_tol_abs = mz_tol[0]

    for i, s in enumerate(relevant_scans):
        intensity, exact_mz = get_max_mz(s, t-mz_tol_abs, t+mz_tol_abs)
        if intensity >= max_i:
            max_i = intensity
            max_idx = i
            max_mz = exact_mz
            max_rt = s.rt_in_seconds
            max_scan_no = s.scan_no
    isos.append((0, t, max_i, max_mz, max_rt, max_scan_no))
    pos = 0
    for t in target_mz[1:]:
        pos += 1
        if pos > max_iso_n:
            break
        max_i = 0
        max_mz = -1
        max_rt = None
        max_scan_no = None
        for scan_idx in range(max_idx-scan_delta, max_idx+scan_delta+1):
            if scan_idx >= 0 and scan_idx < len(relevant_scans):
                s = relevant_scans[scan_idx]
                intensity, exact_mz = get_max_mz(s, t-mz_tol_abs, t+mz_tol_abs)
                if intensity >= max_i:
                    max_i = intensity
                    max_mz = exact_mz
                    max_rt = s.rt_in_seconds
                    max_scan_no = s.scan_no
        isos.append((pos, t, max_i, max_mz, max_rt, max_scan_no))
    return isos


def get_max_mz(scan, mz_min, mz_max):
    sub_peaks = list(
        filter(lambda x: x[0] >= mz_min and x[0] <= mz_max, scan.peaks))
    if len(sub_peaks) == 0:
        return 0.0, -1
    else:
        mz, intensity = zip(*sub_peaks)
        max_i = max(intensity)
        max_mz = np.array(mz)[np.argmax(intensity)]
        return max_i, max_mz


def compute_lipid_kinetics(lipid_name, lipid_dict, file_time_list, mzml_file_objs, parameters):
    rt_range = lipid_dict['rt_range']

    exclude_files = lipid_dict.get('files to exclude', "")
    if len(exclude_files) == 0:
        n_exclude = 0
    else:
        n_exclude = len(exclude_files.split(';'))

    n_iso = lipid_dict.get('max_iso_n', parameters['max_iso_n'])
    mz_tolerance = lipid_dict.get('mz_tolerance', parameters['mz_tolerance'])
    mz_tolerance_units = lipid_dict.get(
        'mz_tolerance_units', parameters['mz_tolerance_units'])
    mz_tol = (mz_tolerance, mz_tolerance_units)

    data_mat = np.zeros((len(file_time_list)-n_exclude, n_iso+1))

    discard_pos = -1
    all_isos = {}
    bespoke_file_time_list = []
    exclude_files = lipid_dict.get('files_to_exclude', "").split(';')
    for fpos, (o, time_val) in enumerate(file_time_list):

        if o in exclude_files:
            print("\t{}, ignoring {}".format(lipid_name, o))
            continue
        bespoke_file_time_list.append((o, time_val))
        mzml_file_obj = mzml_file_objs[o]
        isos = get_iso_intense(mzml_file_obj, rt_range,
                               lipid_dict['formula'],
                               lipid_dict['adduct_type'],
                               #                                mz_tol = (lipid_dict['mz_tolerance_ppm (optional)'],'ppm'),
                               mz_tol=mz_tol,
                               max_iso_n=n_iso,
                               scan_delta=parameters['scan_delta'])

        if fpos == 0:
            # check that intensities are monotonically decreasig. If not, chop.
            for a, (b, _, intensity, _, _, _) in enumerate(isos[:-1]):
                if intensity < isos[a+1][2]:
                    discard_pos = a
                    print(lipid_name, discard_pos)
                    break
        for ipos, (_, _, intensity, _, _, _) in enumerate(isos):
            if discard_pos > -1 and ipos > discard_pos:
                print(lipid_name, fpos, ipos)
                data_mat[fpos, ipos] = 0
            else:
                data_mat[fpos, ipos] = intensity

        all_isos[o] = isos

    if discard_pos > -1:
        data_mat = data_mat[:, :discard_pos+1]

    times = [t for (f, t) in bespoke_file_time_list]

    data_mat /= data_mat.sum(axis=1)[:, None]

    p = fit(times, data_mat, fix_ends=False, make_plot=False)
    k, a0, ai = p

    output_dict = {}
    output_dict['kinetic_parameters'] = (k, a0, ai)
    output_dict['data_mat'] = data_mat
    output_dict['times'] = times
    output_dict['all_isos'] = all_isos

    return output_dict


def create_plot(lipid_name, output_dict, output_filename=None):
    times = output_dict['times']
    data_mat = output_dict['data_mat']
    k, a0, ai = output_dict['kinetic_parameters']

    plt.figure(figsize=(20, 4))
    plt.subplot(1, 3, 1)
    plt.imshow(data_mat, aspect='auto')
    plt.xlabel('isotope')
    plt.ylabel('time')
    plt.yticks(range(len(times)), times)
    plt.title(lipid_name)

    plt.colorbar()
    plt.subplot(1, 3, 2)
    t = np.array(times)
    plt.plot(t, data_mat[:, 0], 'ro')
    plt.plot(t, a0 + (ai-a0)*(1-np.exp(-k*t)))
    plt.title('{}  k: {:.3f}, a0: {:.3f}, ai: {:.3f}'.format(
        lipid_name, k, a0, ai))

    plt_name = '{}_condition_{}_{}.png'.format(
        '_'.join(lipid_name.split()), 1, 'Pos').replace(':', '_')

    all_isos = output_dict['all_isos']
    plt.subplot(1, 3, 3)

    for f, iso_dat in all_isos.items():
        mz = []
        rt = []
        for iso in iso_dat:
            iso_mz = iso[3]
            iso_rt = iso[4]
            if iso_mz > 0:  # it's -1 if not found
                mz.append(iso_mz)
                rt.append(iso_rt)
        plt.scatter(rt, mz, label="{:5.5s}".format(f.split(os.sep)[-1]))
    plt.xlabel('rt')
    plt.ylabel('mz')
    plt.legend()

    if output_filename:
        #         plt.savefig(os.path.join('plots',plt_name))
        #         print("Writing: ",plt_name)
        plt.savefig(output_filename)
        print("Writing: ", output_filename)
    plt.close()


def col2alphabet(col_num):
    prefix = (col_num // 26)
    nextfix = col_num % 26
    if prefix == 0:
        alppre = ''
    else:
        alppre = chr(prefix + 65 - 1)
    alp = chr(nextfix + 65)
    return alppre + alp


def write_xlsx_block(worksheet, list_vals, row_num, start_col):
    for i, v in enumerate(list_vals):
        cell = col2alphabet(start_col+i) + str(row_num)
        worksheet.write(cell, v)


def create_xlsx_output(output_dict, output_filename='test.xlsx'):
    workbook = xlsxwriter.Workbook(
        output_filename, {'nan_inf_to_errors': True})
    for lipid_no, lipid in enumerate(output_dict):
        row_num = 1
        sheet_name = lipid.replace(':', ' ')
        sheet_name = sheet_name.replace('[', '(')
        sheet_name = sheet_name.replace(']', ')')
        worksheet = workbook.add_worksheet(sheet_name)
        cell = col2alphabet(1) + str(row_num)
        worksheet.write(cell, lipid)
        row_num += 1
        write_list = ['Kinetic Parameters', 'k:',
                      output_dict[lipid]['kinetic_parameters'][0],
                      'ai:',
                      output_dict[lipid]['kinetic_parameters'][1],
                      'a0:',
                      output_dict[lipid]['kinetic_parameters'][2]]
        write_xlsx_block(worksheet, write_list, row_num, 1)
        row_num += 2
        write_xlsx_block(worksheet, ['isotope data'], row_num, 1)

        data_mat = output_dict[lipid]['data_mat']
        n_time, n_iso = data_mat.shape

        iso_list = range(n_iso)
        times = output_dict[lipid]['times']
        write_xlsx_block(worksheet, iso_list, row_num, 3)
        row_num += 1
        for i, row in enumerate(data_mat):
            dlist = [times[i]] + list(row)
            write_xlsx_block(worksheet, dlist, row_num, 2)
            row_num += 1

        row_num += 1
        write_xlsx_block(worksheet, ['isotope details', 'file', 'iso number',
                                     'theoretical m/z', 'intensity', 'mz', 'rt', 'scan number'], row_num, 1)
        row_num += 1
        all_isos = output_dict[lipid]['all_isos']
        for filename, iso_data in all_isos.items():
            write_xlsx_block(worksheet, [filename], row_num, 2)
            row_num += 1
            for row in iso_data:
                write_xlsx_block(worksheet, row, row_num, 3)
                row_num += 1

        create_plot(
            lipid, output_dict[lipid], output_filename='temp_{}.png'.format(lipid_no))
        worksheet.insert_image('N4', 'temp_{}.png'.format(lipid_no))
    workbook.close()

    for lipid_no, lipid in enumerate(output_dict):
        os.system('rm temp_{}.png'.format(lipid_no))
