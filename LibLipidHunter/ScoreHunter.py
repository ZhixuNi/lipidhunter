# -*- coding: utf-8 -*-
#
# Copyright (C) 2016-2019  SysMedOs_team @ AG Bioanalytik, University of Leipzig:
# SysMedOs_team: Zhixu Ni, Georgia Angelidou, Mike Lange, Maria Fedorova
# LipidHunter is Dual-licensed
#     For academic and non-commercial use: `GPLv2 License` Please read more information by the following link:
#         [The GNU General Public License version 2] (https://www.gnu.org/licenses/old-licenses/gpl-2.0.en.html)
#     For commercial use:
#         please contact the SysMedOs_team by email.
# Please cite our publication in an appropriate form.
# Ni, Zhixu, Georgia Angelidou, Mike Lange, Ralf Hoffmann, and Maria Fedorova.
# "LipidHunter identifies phospholipids by high-throughput processing of LC-MS and shotgun lipidomics datasets."
# Analytical Chemistry (2017).
# DOI: 10.1021/acs.analchem.7b01126
#
# For more info please contact:
#     Developer Zhixu Ni zhixu.ni@uni-leipzig.de
#     Developer Georgia Angelidou georgia.angelidou@uni-leipzig.de

import os
import re

import pandas as pd

from LibLipidHunter.IsotopeHunter import IsotopeHunter
from LibLipidHunter.AbbrElemCalc import ElemCalc
from LibLipidHunter.LipidNomenclature import NameParserFA
from LibLipidHunter.PanelPlotter import plot_spectra
from LibLipidHunter.PanelPlotter import gen_plot

# TODO: check how the program work for the other lipid classes apart TG
def get_specific_peaks(key_frag_dct, mz_lib, ms2_df, ms2_ppm=100, vendor='waters', exp_mode='LC-MS'):
    ms2_max_i = ms2_df['i'].max()
    ms2_precision = ms2_ppm * 0.000001
    target_frag_df = key_frag_dct['target_frag_df']
    target_nl_df = key_frag_dct['target_nl_df']
    other_frag_df = key_frag_dct['other_frag_df']
    other_nl_df = key_frag_dct['other_nl_df']

    _target_frag_df = pd.DataFrame()
    _target_nl_df = pd.DataFrame()
    _other_frag_df = pd.DataFrame()
    _other_nl_df = pd.DataFrame()

    for _i, _frag_se in target_frag_df.iterrows():

        _frag_mz = _frag_se['EXACTMASS']
        _frag_class = _frag_se['CLASS']
        _frag_label = _frag_se['LABEL']
        _frag_mz_query_code = '%f <= mz <= %f' % (_frag_mz * (1 - ms2_precision), _frag_mz * (1 + ms2_precision))

        _frag_df = ms2_df.query(_frag_mz_query_code)

        if not _frag_df.empty:
            _frag_df = _frag_df.sort_values(by='i', ascending=False)
            _frag_df.loc[:, 'CLASS'] = _frag_class
            _frag_df.loc[:, 'LABEL'] = _frag_label
            _frag_df.loc[:, _frag_label] = 100 * _frag_df['i'] / ms2_max_i
            _target_frag_df = _target_frag_df.append(_frag_df.head(1))

    if not other_frag_df.empty:
        for _i, _frag_se in other_frag_df.iterrows():

            _frag_mz = _frag_se['EXACTMASS']
            _frag_class = _frag_se['CLASS']
            _frag_label = _frag_se['LABEL']
            _frag_mz_query_code = '%f <= mz <= %f' % (_frag_mz * (1 - ms2_precision), _frag_mz * (1 + ms2_precision))

            _frag_df = ms2_df.query(_frag_mz_query_code)

            if not _frag_df.empty:
                _frag_df = _frag_df.sort_values(by='i', ascending=False)
                _frag_df.loc[:, 'CLASS'] = _frag_class
                _frag_df.loc[:, 'LABEL'] = _frag_label
                _frag_df.loc[:, _frag_label] = 100 * _frag_df['i'] / ms2_max_i
                _other_frag_df = _other_frag_df.append(_frag_df.head(1))

    for _i, _nl_se in target_nl_df.iterrows():

        _nl_mz = _nl_se['EXACTMASS']
        _nl_class = _nl_se['CLASS']
        _nl_label = _nl_se['LABEL']
        _nl_mz_query_code = '%f <= mz <= %f' % ((mz_lib - _nl_mz) * (1 - ms2_precision),
                                                (mz_lib - _nl_mz) * (1 + ms2_precision))

        _nl_df = ms2_df.query(_nl_mz_query_code)

        if not _nl_df.empty:
            _nl_df = _nl_df.sort_values(by='i', ascending=False)
            _nl_df.loc[:, 'CLASS'] = _nl_class
            _nl_df.loc[:, 'LABEL'] = _nl_label
            _nl_df.loc[:, _nl_label] = 100 * _nl_df['i'] / ms2_max_i
            _target_nl_df = _target_nl_df.append(_nl_df.head(1))

    for _i, _nl_se in other_nl_df.iterrows():

        _nl_mz = _nl_se['EXACTMASS']
        _nl_class = _nl_se['CLASS']
        _nl_label = _nl_se['LABEL']
        _nl_mz_query_code = '%f <= mz <= %f' % ((mz_lib - _nl_mz) * (1 - ms2_precision),
                                                (mz_lib - _nl_mz) * (1 + ms2_precision))
        _nl_df = ms2_df.query(_nl_mz_query_code)

        if not _nl_df.empty:
            _nl_df = _nl_df.sort_values(by='i', ascending=False)
            _nl_df.loc[:, 'CLASS'] = _nl_class
            _nl_df.loc[:, 'LABEL'] = _nl_label
            _nl_df.loc[:, _nl_label] = 100 * _nl_df['i'] / ms2_max_i
            _other_nl_df = _other_nl_df.append(_nl_df.head(1))

    specific_ion_dct = {}
    if not _target_frag_df.empty:
        specific_ion_dct['TARGET_FRAG'] = _target_frag_df
    if not _target_nl_df.empty:
        specific_ion_dct['TARGET_NL'] = _target_nl_df
    if not _other_frag_df.empty:
        specific_ion_dct['OTHER_FRAG'] = _other_frag_df
    if not _other_nl_df.empty:
        specific_ion_dct['OTHER_NL'] = _other_nl_df

    return specific_ion_dct

def get_all_fa_nl(fa_df, ms2_df, peak_type_lst, lipid_type='LPL'):
    dg_fa_rgx = re.compile(r'\[M-\(*(FA\d{1,2}:\d)')
    # \[M\-(?:(FA\d{1,2}\:\d)|\((FA\d{1,2}\:\d).*\))\+(?:H|Na)\]\+    More strict regular expression
    fa_regexp = re.compile(r'^(FA\d)')
    fa_regexp2 = re.compile(r'(.*\((FA\d)-*.*\).*)')
    m_nl_regexp = re.compile(r'\[[a-zA-Z][\+\-](\((FA\d{1,2}:\d{1,2})[\+\-]*[a-zA-Z0-9]*\))|\[[a-zA-Z][\+\-](([a-zA-Z0-9]*[\-\+]*)*)[\-\+][a-zA-Z0-9]\]')
    name_parser = NameParserFA(lipid_type)
    gen_fragments_dct = name_parser.lipid_fragment_info
    obs_peaks_df = pd.DataFrame()
    bp_i = ms2_df['i'].max()
    # Note: For the PL will not work because of different fragments name.
    # Note: Should put around here the section for the general name of the fragments
    # TODO: (need to figure out this section here)
    # Need to add a control check in case there columns inside the fa-df with nan values (which cause problems)
    for _idx, _fa_se in fa_df.iterrows():
        for peak_typ in peak_type_lst:

            # !! IMPORTANT HERE !!
            # TG use lite_info_df for get_all_fa_nl
            # PL use fa_df for get_all_fa_nl

            # default peak_typ for TG is with FA1/FA2/FA3
            # e.g. ['[M-(FA1)+H]+', '[M-(FA2)+H]+', '[M-(FA3)+H]+']
            # default peak_typ for PL is with NO FA assignment
            # e.g. ['[FA-H]-', '[LPE-H]-', '[LPE-H2O-H]-']

            _q_str = _fa_se[peak_typ]
            _q_tmp_df = ms2_df.query(_q_str).copy()
            # _q_tmp_df.is_copy = False
            peak_typ = peak_typ.replace('_Q', '')
            if not _q_tmp_df.empty:
                _q_tmp_df.loc[:, 'lib_mz'] = _fa_se['%s_MZ' % peak_typ]
                _q_tmp_df.loc[:, 'obs_mz'] = _q_tmp_df['mz']
                _q_tmp_df.loc[:, 'obs_i_r'] = 100 * _q_tmp_df['i'] / bp_i
                _q_tmp_df.loc[:, 'obs_ppm'] = 1e6 * (_q_tmp_df['mz'] - _q_tmp_df['lib_mz']) / _q_tmp_df['lib_mz']
                _q_tmp_df.loc[:, 'obs_ppm'] = _q_tmp_df['obs_ppm'].astype(int)
                _q_tmp_df.loc[:, 'obs_ppm_abs'] = _q_tmp_df['obs_ppm'].abs()
                _q_tmp_df.loc[:, 'obs_abbr'] = _fa_se['%s_ABBR' % peak_typ]
                _q_tmp_df.loc[:, 'obs_type'] = peak_typ
                _q_tmp_df.loc[:, 'obs_label'] = _q_tmp_df['obs_mz']
                _q_tmp_df = _q_tmp_df.round({'obs_mz': 4, 'obs_i_r': 1, 'obs_label': 2})
                _q_tmp_df.loc[:, 'obs_label'] = _q_tmp_df['obs_label'].astype(str)

                # Note: maybe is need to move to another section
                for _k, _v in gen_fragments_dct.items():
                    if peak_typ in _v:
                        _q_tmp_df.loc[:, 'obs_type_g'] = _k
                        break
                # Note: add a control statement


                if re.match(fa_regexp, peak_typ):
                    _fa_abbr_match = re.match(fa_regexp, peak_typ)
                    _fa_abbr_lst = _fa_abbr_match.groups()
                    _q_tmp_df.loc[:, 'fa_abbr'] = _fa_se['%s_ABBR' % _fa_abbr_lst[0]]
                    _q_tmp_df.loc[:, 'fa'] = _fa_abbr_lst[0]
                # elif re.match(m_nl_regexp, peak_typ):
                #     _fa_abbr_match = re.match(m_nl_regexp , peak_typ)
                #     _fa_abbr_lst = _fa_abbr_match.groups()
                #     if _fa_abbr_lst[1] is not None:
                #         _q_tmp_df.loc[:, 'fa_abbr'] = _fa_se['%s_ABBR' % _fa_abbr_lst[1]]
                #     else:
                #         # Note: maybe need some control
                #         _q_tmp_df.loc[:, 'fa_abbr'] = _fa_se['%s_ABBR' % _fa_abbr_lst[2]]
                else:

                    _fa_abbr_match = re.match(m_nl_regexp, _fa_se['%s_ABBR' % peak_typ])
                    _fa_abbr_lst = _fa_abbr_match.groups()
                    _fa_abbr_m2 = re.match(fa_regexp2, peak_typ)
                    if _fa_abbr_lst[1] is not None:
                        _q_tmp_df.loc[:, 'fa_abbr'] = _fa_abbr_lst[1]
                    else:
                        # Note: maybe need some control
                        _q_tmp_df.loc[:, 'fa_abbr'] = _fa_abbr_lst[2]
                    if _fa_abbr_m2:
                        _fa_abbr_lst_g = _fa_abbr_m2.groups()
                        _q_tmp_df.loc[:, 'fa'] = _fa_abbr_lst_g[1]
                    else:
                        # Note: Temporary define it like this
                        _q_tmp_df.loc[:, 'fa'] = 'M'
                obs_peaks_df = obs_peaks_df.append(_q_tmp_df)
            else:
                pass


    if not obs_peaks_df.empty:
        obs_peaks_df.sort_values(by=['obs_abbr', 'i', 'obs_ppm_abs'], ascending=[False, False, True], inplace=True)
        obs_peaks_df.drop_duplicates(subset=['obs_abbr', 'obs_type'], keep='first', inplace=True)
        obs_peaks_df.sort_values(by=['i', 'obs_ppm_abs'], ascending=[False, True], inplace=True)
        obs_peaks_df.reset_index(inplace=True, drop=True)
        obs_peaks_df['obs_rank'] = obs_peaks_df.index + 1
        # return obs_peaks_df.head(10)  # Since all fragments in then should sent more than 10
        return obs_peaks_df
    else:
        # no identification
        return obs_peaks_df



def fragment_rank (w_df, ob_d):

    all_weight_groups_lst = w_df['Group'].unique().tolist()
    obs_df_n = pd.DataFrame()



    for _g in all_weight_groups_lst:
        frag_info_lst = w_df.query('Group == {g}'.format(g=_g)).index.tolist()
        obs_df_s = ob_d.query('obs_type in @frag_info_lst').loc[:, ['i', 'fa_abbr', 'obs_type_g', 'obs_ppm_abs', 'obs_i_r']]
        obs_df_s.drop_duplicates(subset=['fa_abbr', 'obs_type_g'], keep='first', inplace=True)
        obs_df_s.sort_values(by=['i', 'obs_ppm_abs'], ascending=[False, True], inplace=True)
        obs_df_s.reset_index(inplace=True, drop=True)
        obs_df_s['obs_rank'] = obs_df_s.index +1
        obs_df_s['Group'] = _g
        obs_df_s = obs_df_s.set_index(['fa_abbr', 'obs_type_g'])
        obs_df_n = obs_df_n.append(obs_df_s)

    return (obs_df_n)

def prep_rankscore (obs_d, obs_dct, origin_info_df, sliced_info_df, weight_dct, lipid_class='DG'):
    ident_obs_peak_df = pd.DataFrame()
    lipid_c_fa = NameParserFA(lipid_class).lipid_fa_dct[lipid_class]
    print (lipid_c_fa[1])
    _fa_sn_lst = lipid_c_fa[1]


    if "M" in obs_d.loc[:, 'fa'].tolist():
        _fa_sn_lst.append('M')

    for _idx, _lite_se in sliced_info_df.iterrows():
        _lipid_abbr = _lite_se['DISCRETE_ABBR']
        _rank_score = 0.0

        # check if FA contribute to multiple sn --> i needs to be divided accordingly
        _ident_fa_abbr_lst = []
        _ident_peak_count = 0
        _sum_fa_abbr_lst = []
        _sum_fa_abbr_dct = {}
        unique_fa_abbr_dct = {}
        fa_multi_dct = {}  # store multiple FA

        _fa_f = 0
        for _i in _fa_sn_lst:
            print (_i)
            if _i == 'M':
                _fa_abbr = 'H2O'
            else:
                _fa_abbr = _lite_se['{fa}_ABBR'.format(fa=_i)]
                _sum_fa_abbr_lst.append(_fa_abbr)
                if _fa_abbr in _sum_fa_abbr_dct.keys():
                    _sum_fa_abbr_dct[_fa_abbr].append(_i)
                else:
                    _sum_fa_abbr_dct[_fa_abbr] = [_i]
            _obs_df_fa = obs_d.query('fa == @_i and fa_abbr == @_fa_abbr')
            print ('1')
            if not _obs_df_fa.empty:
                # Note: checking if found any type of fragments, later will included only specific ones
                _fa_f = _fa_f + 1
                print('2')
                _obs_type_frag_lst = _obs_df_fa['obs_type'].tolist()
                print('3')
                for _idx2,_f_s in _obs_df_fa.iterrows():
                    _f_w = _f_s['obs_type']
                    print (_f_w)
                    if _f_w in weight_dct.index.tolist():
                        origin_info_df.at[_idx, '{f}_i'.format(f=_f_w)] = _f_s['i']
                        origin_info_df.at[_idx, '{f}_i_per'.format(f=_f_w)] = _f_s['obs_i_r']
                        # Note: maybe in this section we should put kind of control if the fragments are in the df
                        obs_r = obs_dct.loc[(_fa_abbr, _f_s['obs_type_g']), 'obs_rank']
                        origin_info_df.at[_idx, '{f}_RANK'.format(f=_f_w)] = obs_r
                        origin_info_df.at[_idx, '{f}_SCORE'.format(f=_f_w)] = ((11 - obs_r) * 0.1 * weight_dct.loc[_f_w, 'Weight'])
                        if _i == 'M':
                            _rank_score += origin_info_df.loc[_idx, '{f}_SCORE'.format(f=_f_w)]
                            _ident_peak_count += 1
        print ('come on')
        origin_info_df.at[_idx, 'RESI_CHK'] = _fa_f
        origin_info_df.at[_idx, 'total_fa'] = len(_fa_sn_lst)

        _unique_fa_abbr_lst = list(set(_sum_fa_abbr_lst))
        print ('well where is the issue')
        for _fa_abbr in _unique_fa_abbr_lst:
            print (_fa_abbr)
            _fa_count = _sum_fa_abbr_lst.count(_fa_abbr)
            _obs_df_fa_2 = obs_d.copy()
            _fa_lst = _sum_fa_abbr_dct[_fa_abbr]

            if _fa_count > 1:

                _obs_df_fa_2.i.loc[(_obs_df_fa_2['fa_abbr'] == _fa_abbr)] = _obs_df_fa_2.i.loc[(_obs_df_fa_2['fa_abbr'] == _fa_abbr)] / _fa_count
                _obs_df_fa_2.obs_i_r.loc[(_obs_df_fa_2['fa_abbr'] == _fa_abbr)] = _obs_df_fa_2.obs_i_r.loc[(
                            _obs_df_fa_2['fa_abbr'] == _fa_abbr)] / _fa_count

                # TODO: figure out what is going when this going inside the fragment_rank
                _obs_df_fa_r = fragment_rank(weight_dct, _obs_df_fa_2)
                # Issue: Does not replace correctly the values which it should
                for _fa_n in _fa_lst:
                    _obs_df_fa = _obs_df_fa_2.query('fa == @_fa_n & fa_abbr == @_fa_abbr')
                    if not _obs_df_fa.empty:
                        # Note: checking if found any type of fragments, later will included only specific ones
                        for _idx2, _f_s in _obs_df_fa.iterrows():
                            _f_w = _f_s['obs_type']
                            if _f_w in weight_dct.index.tolist():
                                origin_info_df.loc[_idx, '{f}_i'.format(f=_f_w)] = _obs_df_fa_r.loc[(_fa_abbr, _f_s['obs_type_g']), 'i']
                                origin_info_df.loc[_idx, '{f}_i_per'.format(f=_f_w)] = _obs_df_fa_r.loc[(_fa_abbr, _f_s['obs_type_g']),'obs_i_r']
                                origin_info_df.loc[_idx, '{f}_RANK'.format(f=_f_w)] = _obs_df_fa_r.loc[(_fa_abbr, _f_s['obs_type_g']),'obs_rank']
                                origin_info_df.loc[_idx, '{f}_SCORE'.format(f=_f_w)] = (
                                            (11 - _obs_df_fa_r.loc[(_fa_abbr, _f_s['obs_type_g']),'obs_rank']) * 0.1 * weight_dct.loc[_f_w, 'Weight'])
                                _rank_score += origin_info_df.loc[_idx, '{f}_SCORE'.format(f=_f_w)]
                                _ident_peak_count += 1
            else:

                _obs_df_fa3 = _obs_df_fa_2.query('fa == @_fa_lst[0] & fa_abbr == @_fa_abbr')

                for _f_s3 in _obs_df_fa3['obs_type'].values.tolist():
                    _rank_score += origin_info_df.loc[_idx, '{f}_SCORE'.format(f=_f_s3)]

                _ident_peak_count += 1


                # _obs_df_fa

        # unique_fa_abbr_dct[_fa_abbr] = sorted(_fa_abbr_sn_lst)
        print ('Are you out of your mind')
        _ident_fa_abbr_lst = list(set(_ident_fa_abbr_lst))
        # print(_lipid_abbr, 'rank_score', _rank_score)
        # Question: Why to we need two control check what is the reason for this?
        # if _rank_score > 0 and len(_ident_fa_abbr_lst) > 0:
        if _rank_score > 0 :
            # print(_lipid_abbr, _rank_score)
            origin_info_df.at[_idx, 'RANK_SCORE'] = _rank_score
            origin_info_df.at[_idx, 'OBS_PEAKS'] = int(_ident_peak_count)
            origin_info_df.at[_idx, 'OBS_RESIDUES'] = len(_sum_fa_abbr_lst)
        # else:
        #     ident_obs_peak_df = pd.DataFrame({'i': [], 'mz': [], 'lib_mz': [], 'obs_mz': [], 'obs_i_r': [],
        #                             'obs_ppm': [], 'obs_ppm_abs': [], 'obs_abbr': [], 'obs_label': [],
        #                             'fa_abbr': [], 'obs_rank': [], 'obs_type_calc': []})
    # Note: Currently ident_obs_peak_df is empty. Untill will now why it is needed and to avoid a few of the possible errors which may occur
    return ident_obs_peak_df, origin_info_df

def calc_rankscore(obs_df, lite_info_df, lipid_class, weight_dct, rankscore_filter, all_sn=True):
    # print('Start to calc for rank score...')

    # Take by one each observation and then check for the the different position of the observed fragment
    # ident_peak_dct = {'discrete_abbr': [], 'obs_label': [], 'i': [], 'mz': [], 'obs_abbr': [], 'obs_rank_type': [],
    #                   'obs_rank': []}
    lite_info_df2 = pd.DataFrame()
    lite_info_df2['DISCRETE_ABBR'] = lite_info_df['DISCRETE_ABBR']
    lite_info_df2['MS2_PR_mz'] = lite_info_df['MS2_PR_mz']
    lite_info_df2['DDA_rank'] = lite_info_df['DDA_rank']
    lite_info_df2['MS1_XIC_mz'] = lite_info_df['MS1_XIC_mz']
    lite_info_df2['ppm'] = lite_info_df['ppm']
    lite_info_df2['MS1_obs_mz'] = lite_info_df['MS1_obs_mz']
    lite_info_df2['Lib_mz'] = lite_info_df['Lib_mz']
    lite_info_df2['BULK_ABBR'] = lite_info_df['BULK_ABBR']
    lite_info_df2['scan_time'] = lite_info_df['scan_time']
    lite_info_df2['RANK_SCORE'] = 0.0
    lite_info_df2['OBS_RESIDUES'] = 0

    if lipid_class in ['TG']:
        resi_site_lst = ['FA1_ABBR', 'FA2_ABBR', 'FA3_ABBR']
    elif lipid_class in ['PA', 'PC', 'PE', 'PG', 'PS', 'PI', 'PIP', 'DG', 'SM']:
        resi_site_lst = ['FA1_ABBR', 'FA2_ABBR']
    elif lipid_class in ['LPA', 'LPC', 'LPE', 'LPG', 'LPS', 'LPI', 'LPIP']:
        resi_site_lst = ['FA1_ABBR']
    elif lipid_class in ['MG']:
        resi_site_lst = ['FA1_ABBR']
    else:
        resi_site_lst = ['FA1_ABBR', 'FA2_ABBR']

    obs_fa_lst = []
    obs_df_n = fragment_rank(weight_dct, obs_df)

    if not lite_info_df.empty:
        # Comment: for know give the lite_info_df untill figure out how to handle and replasce the above code
        ident_peak_df, lite_info_df2 = prep_rankscore(obs_df, obs_df_n, lite_info_df2, lite_info_df, weight_dct,
                                                      lipid_class=lipid_class)
    else:
        ident_peak_df = pd.DataFrame()

    # Question: how the case where there O- and P- avoiding the error since in those case we do not see any fa fragment?
    # Note: Temporary solution
    if lipid_class != 'Cer':
        lite_info_df2['RANK_SCORE'] = (lite_info_df2['RANK_SCORE'] * lite_info_df2['OBS_RESIDUES']/ lite_info_df2['total_fa']).round(2)

    if not lite_info_df2.empty:
        lite_info_df2 = lite_info_df2[lite_info_df2['RANK_SCORE'] >= rankscore_filter]
    # Original code: will activate later when there will understand what for
    # if not ident_peak_df.empty:
    #     ident_peak_df = ident_peak_df[ident_peak_df['discrete_abbr'].isin(lite_info_df2['DISCRETE_ABBR']
    #                                                                       .values.tolist())]
    # if not lite_info_df.empty:
    #     print(lite_info_df[['BULK_ABBR', 'DISCRETE_ABBR', 'RANK_SCORE', 'scan_time', 'OBS_PEAKS']])

    return ident_peak_df, lite_info_df2
#

def get_rankscore(fa_df, master_info_df, abbr_bulk, charge, ms2_df, _ms2_idx, lipid_class, weight_dct, core_count,
                  rankscore_filter=27.5, all_sn=True):
    frag_reg_expr = re.compile(r'^(.*\[.*\][-,\+])_Q$')
    # Question: Why we do not move this query before entering the function (georgia: 14.2.2019)

    lite_info_df = master_info_df.query('BULK_ABBR == "%s" and spec_index == %f' % (abbr_bulk, _ms2_idx)).copy()

    # lite_info_df.is_copy = False
    lite_info_df['RANK_SCORE'] = 0
    lite_info_col_lst = list(lite_info_df.columns.values)
    lite_groups = list(filter(frag_reg_expr.match, lite_info_col_lst))
    obs_fa_frag_df_g = pd.DataFrame(get_all_fa_nl(lite_info_df, ms2_df, lite_groups, lipid_class))

    obs_dct = {}
    frag_lst_dg = []
    frag_lst_dg_w = []

    if obs_fa_frag_df_g.shape[0] > 0:
        post_ident_peak_df, lite_info_df = calc_rankscore(obs_fa_frag_df_g, lite_info_df, lipid_class,
                                                          weight_dct, rankscore_filter, all_sn=True)
    else:
        post_ident_peak_df = pd.DataFrame()
        lite_info_df = pd.DataFrame()

    matched_checker = 0

    if isinstance(lite_info_df, pd.DataFrame):
        if not lite_info_df.empty:
            lite_info_df.sort_values(by=['RANK_SCORE', 'DISCRETE_ABBR'], ascending=[False, True], inplace=True)
            lite_info_df.reset_index(drop=True, inplace=True)
            matched_checker = 1
            # Original code: Unknown the reason for this section
            # if not post_ident_peak_df.empty:
            #     matched_checker = 1
            #     checked_abbr_lst = lite_info_df['DISCRETE_ABBR'].values.tolist()
            #     post_ident_peak_df = post_ident_peak_df[post_ident_peak_df['discrete_abbr'].isin(checked_abbr_lst)]
            #     post_ident_peak_df.sort_values(by='mz', inplace=True)
            #     post_ident_peak_df.reset_index(drop=True, inplace=True)

        # TODO: need to change the abbreviation OBS_FA and OBS_LYSO to more general one
        # Note: OBS_FA indicates the information going to the fith graph in the figure
    out_form_dict = NameParserFA(lipid_class).lipid_output_format[lipid_class][charge]
    obs_info_dct = {'INFO': lite_info_df, 'OBS_FA': out_form_dict['OBS_FA'], 'OBS_LYSO': out_form_dict['OBS_LYSO'],
                        'IDENT': obs_fa_frag_df_g}


    # TODO: check the output from the original program and think how to proceed from here (21.10.2019)
    return matched_checker, obs_info_dct


def ions_control_check(ch_mode, abbr_bulk, ms1_pr, _ms1_df, core_c ):
    # ion flag will be use as an indicator if our precursor correspond to the proposed ion adduct or not
    _mz_amm_iso_flag = ''
    _mz_amm_iso_flag2 = ''
    isotope_hunter = IsotopeHunter()
    get_element = ElemCalc()

    if ch_mode in ['[M+H]+', '[M+Na]+']:
        _mz_amm_formula, _mz_amm_elem_dct = get_element.get_formula(abbr_bulk, charge='neutral')
        _mz_amm_elem_dct['N'] += 1
        _mz_amm_elem_dct['H'] += 4
        _mz_df_amm = pd.DataFrame()
        _mz_amm_mz = get_element.get_exactmass(_mz_amm_elem_dct)
        _mz_amm_mz2, _mz_amm_form, _mz_amm_Na_mz2, _mz_amm_Na_form = get_element.get_NH3_pos_mode(
            ch_mode, ms1_pr, _mz_amm_elem_dct)
        _frag_mz_query_code = '%f <= mz <= %f' % (_mz_amm_mz2 - 0.2, _mz_amm_mz2 + 0.2)
        # TODO (georgia.angelidou@uni-leipzig.de): Need to check the reason why we do not get any output
        _frag_mz_query_code2 = '%f <= mz <= %f' % (_mz_amm_Na_mz2 - 0.2, _mz_amm_Na_mz2 + 0.2)

        if not _ms1_df.query(_frag_mz_query_code).empty:
            # print('Go and check line 613 from ScoreHunter.py to see what is going on')
            _mz_df_amm = _ms1_df.query(_frag_mz_query_code)
            _mz_df_amm.reset_index(inplace=True, drop=True)
            _mz_amm_i = _mz_df_amm.loc[0, 'i']
            # if charge_mode in ['[M+H]+']:
            amm_formula = 'C' + str(_mz_amm_elem_dct['C']) + 'H' + str(
                _mz_amm_elem_dct['H']) + 'O' + str(_mz_amm_elem_dct['O']) + 'N'
            _mz_amm_iso_flag2 = isotope_hunter.get_isotope_fragments(_mz_amm_mz2, _mz_amm_i,
                                                                     _mz_amm_form, _ms1_df,
                                                                     core_c)
            # else:
            #     _mz_amm_iso_flag2 = 0
        else:
            _mz_amm_i = 0
            _mz_amm_iso_flag2 = 1
        if not _ms1_df.query(_frag_mz_query_code2).empty:
            _mz_df_amm_Na = _ms1_df.query(_frag_mz_query_code2)
            _mz_df_amm_Na.reset_index(inplace=True, drop=True)
            _mz_amm_Na_i = _mz_df_amm_Na.loc[0, 'i']
            # if charge_mode in ['[M+Na]+']:
            #     amm_formula = 'C' + str(_mz_amm_elem_dct['C']) + 'H' + str(
            #         _mz_amm_elem_dct['H']) + 'O' + str(_mz_amm_elem_dct['O']) + 'N'
            _mz_amm_iso_Na_flag2 = isotope_hunter.get_isotope_fragments(_mz_amm_Na_mz2,
                                                                        _mz_amm_Na_i,
                                                                        _mz_amm_Na_form,
                                                                        _ms1_df,
                                                                        core_c)
            # else:
            #     _mz_amm_iso_Na_flag2 = 0
        else:
            _mz_amm_Na_i = 0
            _mz_amm_iso_Na_flag2 = 0

        if _mz_amm_Na_i > _mz_amm_i and ch_mode in ['[M+H]+'] \
                and _mz_amm_iso_flag2 == 0 and _mz_amm_iso_Na_flag2 != 1:
            _mz_amm_iso_flag = 1
        elif _mz_amm_i > _mz_amm_Na_i and ch_mode in ['[M+Na]+'] \
                and _mz_amm_iso_Na_flag2 == 0 and _mz_amm_iso_flag2 != 1:
            _mz_amm_iso_flag = 1
        elif (_mz_amm_iso_flag2 == 1 and ch_mode in ['[M+H]+']) or (
                _mz_amm_iso_Na_flag2 == 1 and ch_mode in ['[M+Na]+']):
            _mz_amm_iso_flag = 1
        # elif _mz_amm_Na_i > 0 and _mz_amm_i == 0 and _mz_amm_iso_flag2 != 1:
        #     _mz_amm_iso_flag =1
        else:
            _mz_amm_iso_flag = 0

    else:
        _mz_amm_iso_flag = 0
    return(_mz_amm_iso_flag)


def get_lipid_info(param_dct, fa_df, checked_info_df, checked_info_groups, core_list, usr_weight_df,
                   key_frag_dct, core_spec_dct, xic_dct, core_count, save_fig=True, os_type='windows', queue=None):
    core_count = 'Core_#%i' % core_count

    usr_lipid_type = param_dct['lipid_class']
    charge_mode = param_dct['charge_mode']
    output_folder = param_dct['img_output_folder_str']
    usr_ms2_th = param_dct['ms2_th']
    usr_ms1_precision = param_dct['ms_ppm'] * 1e-6
    usr_ms2_ppm = param_dct['ms2_ppm']
    usr_ms2_info_th_p = param_dct['ms2_infopeak_threshold']

    usr_isotope_score_filter = param_dct['isotope_score_filter']
    usr_rankscore_filter = param_dct['rank_score_filter']

    img_typ = param_dct['img_type']
    img_dpi = param_dct['img_dpi']
    usr_vendor = param_dct['vendor']
    exp_mode = param_dct['experiment_mode']
    usr_fast_isotope = param_dct['fast_isotope']
    usr_tag_all_sn = param_dct['tag_all_sn']

    hunter_start_time_str = param_dct['hunter_start_time']
    isotope_hunter = IsotopeHunter()

    # score_calc = ScoreGenerator(param_dct, usr_weight_df, usr_key_frag_df, usr_lipid_type,
    #                             checked_info_df, ion_charge=charge_mode, ms2_ppm=usr_ms2_ppm)

    # usr_weight_dct = usr_weight_df.to_dict(orient='index')

    tmp_df = pd.DataFrame()

    img_plt_lst = []

    for group_key in core_list:
        _subgroup_df = checked_info_groups.get_group(group_key)
        # TODO (georgia.angelidou@uni-leipzig.de): here a control is need to be done since maybe there some problems If we have one less DB and one more C maybe can cause some problems or maybe not
        # Keep it in mind
        _usr_abbr_bulk_lst = list(set(_subgroup_df['BULK_ABBR'].values.tolist()))
        # TODO (georgia.angelidou@uni-leipzig.de): Here should be the control for the ms values to avoid problems
        usr_spec_info_dct = core_spec_dct[group_key]
        # Note: in the current table there is also another coloumn with the formula info (Formula) (georgia: 14.2.2019)
        key_type = '%s_FORMULA' % _subgroup_df['Ion'].iloc[0]
        # Question: why do we get only the first entry?
        # Example: row1 TG(16:1_16:1_16:2) and row2 TG(16:0_16:2_16:2) what habbens with the second one
        _samemz_se = _subgroup_df.loc[:, _subgroup_df.columns.isin(
            list(('scan_time', 'DDA_rank', 'scan_number', 'Lib_mz', 'FORMULA', 'Ion', key_type, 'Lib_mz', 'MS1_XIC_mz',
                  'MS2_PR_mz')))].iloc[0]
        _usr_ms2_rt = _samemz_se['scan_time']

        _usr_ms2_dda_rank = _samemz_se['DDA_rank']
        _usr_ms2_scan_id = _samemz_se['scan_number']
        _usr_mz_lib = _samemz_se['Lib_mz']
        _usr_formula = _samemz_se['FORMULA']
        _usr_charge = _samemz_se['Ion']
        # DONE: replace '%s_FORMULA' % _usr_charge with key_type (georgia: 14.2.2019)
        # _usr_formula_charged = _samemz_se['%s_FORMULA' % _usr_charge]
        _usr_formula_charged = _samemz_se[key_type]
        _usr_ms2_pr_mz = _samemz_se['Lib_mz']
        _obs_ms2_pr_mz = _samemz_se['MS2_PR_mz']
        _ms1_pr_i = usr_spec_info_dct['ms1_i']
        _ms1_pr_mz = usr_spec_info_dct['ms1_mz']
        _ms1_df = usr_spec_info_dct['ms1_df']
        _ms2_df = usr_spec_info_dct['ms2_df']
        _ms2_idx = usr_spec_info_dct['_ms2_spec_idx']
        _ms1_rt = usr_spec_info_dct['ms1_rt']

        print(core_count, '[INFO] --> ', _usr_ms2_rt, _ms1_pr_mz, _usr_formula_charged)
        # _mz_amm_flag  = isotope_hunter.get_isotope_fragments(_ms1_df, )

        # use the max threshold from abs & relative intensity settings
        if 'i' in _ms2_df.columns.values.tolist():
            _ms2_max_i = _ms2_df['i'].max()
            ms2_threshold = max(usr_ms2_th, _ms2_max_i * usr_ms2_info_th_p)
            _score_ms2_df = _ms2_df.query('i > %f' % ms2_threshold)
        else:
            _ms2_df = pd.DataFrame()
            _score_ms2_df = pd.DataFrame()

        scan_time_chk = False
        # Control to make sure that ms1 is before ms2
        if _ms1_rt > _usr_ms2_rt:
            print(core_count,
                  '[WARNING] !!! MS and MS/MS mismatch MS1_RT {ms1_rt}> MS2_RT {ms2_rt}...'
                  .format(ms1_rt=_ms1_rt, ms2_rt=_usr_ms2_rt))
        else:
            scan_time_chk = True

        spec_df_chk = False
        # Control check if any ms1 or ms2 is empty or not
        if not _ms1_df.empty and not _ms2_df.empty:
            spec_df_chk = True
        elif _ms1_df.empty:
            print(core_count, '[WARNING] !!! MS1 spectrum is empty...')
        elif _ms2_df.empty:
            print(core_count, '[WARNING] !!! MS/MS spectrum is empty...')
        else:
            pass

        if _ms1_pr_mz > 0.0 and _ms1_pr_i > 0.0 and spec_df_chk is True and scan_time_chk is True:
            print(core_count, '[INFO] --> Best PR on MS1: %f' % _ms1_pr_mz)

            isotope_score_info_dct = isotope_hunter.get_isotope_score(_ms1_pr_mz, _ms1_pr_i,
                                                                      _usr_formula_charged, _ms1_df, core_count,
                                                                      ms1_precision=usr_ms1_precision,
                                                                      isotope_number=2,
                                                                      only_c=usr_fast_isotope,
                                                                      score_filter=usr_isotope_score_filter)

            isotope_score = isotope_score_info_dct['isotope_score']
            _mz_amm_iso_flag = ""
            print(core_count, '[SCORE] === isotope_score: %f' % isotope_score)
            if isotope_score >= usr_isotope_score_filter:
                print(core_count, '[INFO] --> isotope_check PASSED! >>> >>> >>>')
                # print(core_count, '>>> >>> >>> >>> Entry Info >>> >>> >>> >>> ')
                _samemz_se.at['MS1_obs_mz'] = _ms1_pr_mz
                _exact_ppm = 1e6 * (_ms1_pr_mz - _usr_mz_lib) / _usr_mz_lib
                _samemz_se.at['ppm'] = _exact_ppm
                _samemz_se.at['abs_ppm'] = abs(_exact_ppm)
                print(core_count, '[INFO] --> Proposed_bulk_structure can be:', _usr_abbr_bulk_lst)
                for _usr_abbr_bulk in _usr_abbr_bulk_lst:
                    # TODO (georgia.angelidou@uni-leipzig.de): Here is where the check for the ammonium adduct should be inlude before the rank score
                    print(core_count, '[INFO] --> Now check_proposed_structure:', _usr_abbr_bulk)
                    # Note: maybe can create a seperate function for the TG.
                    # Note: 2 also ceramide in [M+H]+ this part will be not valid
                    if charge_mode in ['[M+H]+', '[M+Na]+'] and usr_lipid_type == 'TG':
                        _mz_amm_iso_flag = ions_control_check(charge_mode, _usr_abbr_bulk, _ms1_pr_mz, _ms1_df, core_count)
                    else:
                        _mz_amm_iso_flag = 0


                    if _mz_amm_iso_flag == 0:
                        matched_checker, obs_info_dct = get_rankscore(fa_df, checked_info_df, _usr_abbr_bulk,
                                                                      charge_mode,
                                                                      _score_ms2_df, _ms2_idx, usr_lipid_type,
                                                                      usr_weight_df, core_count,
                                                                      rankscore_filter=usr_rankscore_filter,
                                                                      all_sn=usr_tag_all_sn)
                        if matched_checker > 0:
                            # Note: we can replace the obs_info_df and not replace it and use the original abbreviation (georgia: 14.2.2019)
                            obs_info_df = obs_info_dct['INFO']
                            print (obs_info_df.iloc[0])
                            rank_score = obs_info_df['RANK_SCORE'].values.tolist()
                        else:
                            obs_info_dct = {}
                            rank_score = []
                            obs_info_df = pd.DataFrame()
                    else:
                        matched_checker = 0
                        obs_info_dct = {}
                        rank_score = []
                        obs_info_df = pd.DataFrame()

                    if matched_checker > 0:
                        if len(key_frag_dct) > 0:
                            specific_dct = get_specific_peaks(key_frag_dct, _usr_mz_lib, _score_ms2_df,
                                                              ms2_ppm=usr_ms2_ppm, vendor=usr_vendor,
                                                              exp_mode=exp_mode)
                        else:
                            specific_dct = {}

                        specific_ion_count = 0
                        unspecific_ion_count = 0

                        if 'TARGET_FRAG' in list(specific_dct.keys()):
                            specific_ion_count += specific_dct['TARGET_FRAG'].shape[0]
                        if 'TARGET_NL' in list(specific_dct.keys()):
                            specific_ion_count += specific_dct['TARGET_NL'].shape[0]
                        if 'OTHER_FRAG' in list(specific_dct.keys()):
                            unspecific_ion_count += specific_dct['OTHER_FRAG'].shape[0]
                        if 'OTHER_NL' in list(specific_dct.keys()):
                            unspecific_ion_count += specific_dct['OTHER_NL'].shape[0]

                        print(core_count, '[SCORE] === Rank_score: --> passed', rank_score)
                        print(core_count, '\n',
                              obs_info_df[['BULK_ABBR', 'DISCRETE_ABBR', 'RANK_SCORE', 'scan_time', 'OBS_RESIDUES']]
                              )

                        # format abbr. for file names
                        _save_abbr_bulk = _usr_abbr_bulk
                        _save_abbr_bulk = _save_abbr_bulk.replace(r'(', r'[')
                        _save_abbr_bulk = _save_abbr_bulk.replace(r')', r']')
                        _save_abbr_bulk = _save_abbr_bulk.replace(r'<', r'[')
                        _save_abbr_bulk = _save_abbr_bulk.replace(r'>', r']')
                        _save_abbr_bulk = _save_abbr_bulk.replace(r':', r'-')
                        _save_abbr_bulk = _save_abbr_bulk.replace(r'@', r'-')
                        _save_abbr_bulk = _save_abbr_bulk.replace('\\', r'_')
                        _save_abbr_bulk = _save_abbr_bulk.replace(r'/', r'_')

                        img_name_core = ('/%.4f_rt%.3f_DDAtop%.0f_scan%.0f_%s.%s'
                                         % (_usr_ms2_pr_mz, _usr_ms2_rt, _usr_ms2_dda_rank,
                                            _usr_ms2_scan_id, _save_abbr_bulk, img_typ)
                                         )
                        img_name = (output_folder +
                                    r'/LipidHunter_Results_Figures_%s'
                                    % hunter_start_time_str + img_name_core)

                        # print(core_count, '==> check for output -->')
                        # Question: why we keep all of this without changing them. Consider to remove them and reduce the amount of information
                        # TODO: remove all the unnecessary information from the list
                        obs_info_df['Proposed_structures'] = _usr_abbr_bulk
                        obs_info_df['Bulk_identification'] = _usr_abbr_bulk
                        # obs_info_df['Discrete_identification'] = _usr_abbr_bulk
                        obs_info_df['Formula_neutral'] = _usr_formula
                        obs_info_df['Formula_ion'] = _usr_formula_charged
                        obs_info_df['Charge'] = _usr_charge
                        obs_info_df['MS1_obs_mz'] = _ms1_pr_mz
                        obs_info_df['MS1_obs_i'] = '%.2e' % float(_ms1_pr_i)
                        # Note: it is already present in the dataframe (georgia: 15.2.2019)
                        obs_info_df['Lib_mz'] = _usr_mz_lib
                        obs_info_df['MS2_scan_time'] = _usr_ms2_rt
                        obs_info_df['DDA#'] = _usr_ms2_dda_rank
                        # Note: it is already present in the dataframe (georgia: 15.2.2019)
                        obs_info_df['MS2_PR_mz'] = _obs_ms2_pr_mz
                        obs_info_df['Scan#'] = _usr_ms2_scan_id
                        obs_info_df['ISOTOPE_SCORE'] = isotope_score
                        obs_info_df['#Specific_peaks'] = specific_ion_count
                        obs_info_df['#Unspecific_peaks'] = unspecific_ion_count
                        # Note: there is already a column with this name (georgia: 15.2.2019)
                        obs_info_df['ppm'] = _exact_ppm
                        obs_info_df['img_name'] = img_name_core[1:]

                        usr_spec_info_dct['ms1_mz'] = isotope_score_info_dct['obs_pr_mz']
                        usr_spec_info_dct['ms1_i'] = isotope_score_info_dct['obs_pr_i']

                        tmp_df = tmp_df.append(obs_info_df)
                        print (obs_info_df.iloc[0])
                        if save_fig is True:
                            ms1_xic_mz = _samemz_se['MS1_XIC_mz']
                            img_param_dct = {
                                'abbr': _usr_abbr_bulk,
                                'mz_se': _samemz_se,
                                'xic_df': xic_dct[ms1_xic_mz],
                                'ident_info_dct': obs_info_dct,
                                'spec_info_dct': usr_spec_info_dct,
                                'isotope_score_info_dct': isotope_score_info_dct,
                                'specific_dct': specific_dct,
                                'formula_charged': _usr_formula_charged,
                                'charge': _usr_charge,
                                'save_img_as': img_name
                            }

                            img_plt_lst.append(img_param_dct.copy())

                        else:
                            pass
        else:
            pass
    if not tmp_df.empty:
        print(core_count, '[INFO] --> Size of the identified LPP_df %i, %i' % (tmp_df.shape[0], tmp_df.shape[1]))
        tmp_df.reset_index(drop=True, inplace=True)
        tmp_df.index += 1

    else:
        print(core_count, '[WARNING] !!! Size of the identified LPP_df == 0')
        tmp_df = 'error'
    if save_fig is True:
        # print('img_plt_lst', len(img_plt_lst))
        pass
    else:
        img_plt_lst = []

    r_lst = (tmp_df, img_plt_lst)

    if os_type == 'linux_multi':
        queue.put(r_lst)
    else:
        return r_lst
