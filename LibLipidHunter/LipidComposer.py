# -*- coding: utf-8 -*-
#
# Copyright (C) 2016-2017  SysMedOs_team @ AG Bioanalytik, University of Leipzig:
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
#     SysMedOs_team: oxlpp@bbz.uni-leipzig.de
#     Developer Zhixu Ni zhixu.ni@uni-leipzig.de
#     Developer Georgia Angelidou georgia.angelidou@uni-leipzig.de

from __future__ import division
from __future__ import print_function

import itertools

import pandas as pd
import numpy as np
import re
from natsort import natsorted, ns

try:
    from LibLipidHunter.LipidNomenclature import NameParserFA
    from LibLipidHunter.AbbrElemCalc import ElemCalc
    from LibLipidHunter.ParallelFunc import ppm_window_para
except ImportError:  # for python 2.7.14
    from LipidNomenclature import NameParserFA
    from AbbrElemCalc import ElemCalc
    from ParallelFunc import ppm_window_para


class LipidComposer:

    def __init__(self):
        # Question: For what this is used. Remember to check it out (georgia: 23.1.2019)
        self.lipid_hg_lst = ['PA', 'PC', 'PE', 'PG', 'PI', 'PS', 'PIP', 'TG']

        self.lipid_hg_elem_dct = ElemCalc().lipid_hg_elem_dct

        # Note: Already define in the AbbrElemCalc. Will be replace (georgia: 14.2.2019)
        self.glycerol_bone_elem_dct = ElemCalc().glycerol_bone_elem_dct
        self.link_o_elem_dct = {'O': -1, 'H': 2}
        self.link_p_elem_dct = {'O': -1}

        self.elem_dct = NameParserFA().elem_dct

        # Note: need for 2 different functions (georgia: 23.1.2019)
        self.all_lipid_class_list = ['PL', 'PA', 'PC', 'PE', 'PG', 'PI', 'PS', 'SM', 'LPL', 'LPA', 'LPC', 'LPE', 'LPG',
                                'LPI', 'LPS', 'TG', 'DG', 'Cer']



    def calc_fa_df2(self, lipid_class, fa_df):
        # This function creates a list for all the different FA for the different position
        # Will be use in a later step to get all the possible structures that can be combine
        sn_units_lst = []
        # Note: the calc_fa_df was compine in this function (georgia: 15.2.2019)
        header_lst = fa_df.columns.values.tolist()
        lipid_class_dct = NameParserFA(lipid_class).lipid_fa_dct[lipid_class]
        #print (NameParserFA(lipid_class).lipid_fa_dct[lipid_class])
        if lipid_class_dct[0] in header_lst and 'FATTYACID' in header_lst:

            for k in lipid_class_dct[1]:
                q_str = '{cl} == "T" and {fa} == "T"'.format(cl=lipid_class_dct[0], fa=k)
                #print(q_str)
                fa_lst = fa_df.query(q_str)['FATTYACID'].tolist()
                spb_str = re.compile(r'{f}'.format(f=lipid_class_dct[3]))
                spb_lst = list(filter(spb_str.match, fa_lst))
                # Note: this it will work as long the sphingoid base is define in the FA1 and the FA for sphingoid lipids in the FA2
                # Note2: if the user do a mistake then there will be probles. Consider different filters.
                # Note3: Posible filter to check the sphingoid base
                if len(spb_lst) >= 1:
                    sn_units_lst.append(spb_lst)

                else:
                    fa_str = 'FA'
                    fa_lst = list(filter(lambda x: fa_str in x, fa_lst))
                    sn_units_lst.append(fa_lst)
            # # Note: it can be skip and do this in a later step (georgia: 18.2.2019)
            # if NameParserFA().lipid_fa_dct[lipid_class][0] == 'LPL':
            #     sn_units_lst.append('')

        return sn_units_lst

    def calc_fa_query(self, lipid_class, ion_mode, fa_whitelist, ms2_ppm=100):

        usr_fa_df = pd.read_excel(fa_whitelist)
        usr_fa_df = usr_fa_df.fillna(value='F')
        tmp_columns = usr_fa_df.columns.tolist()

        usr_fa_df.columns = usr_fa_df.columns.str.upper()

        # COMPINE (georgia):  the below section was rearrange to one 23.01.19
        #  Note: is more like a control check (georgia:23.1.2019)
        if lipid_class in self.all_lipid_class_list:
            if lipid_class in tmp_columns:
                pass
            elif 'PL' in tmp_columns or 'LPL' in tmp_columns:
                pass
            else:
                return False
        else:
            return False


        sn_units_lst = self.calc_fa_df2(lipid_class, usr_fa_df)  # Return a list with list of the FA for each sn position
        fa_abbr_lst = []
        # For PL lem(sn_units_lst) = 2 and for TG len(sn_units_lst) = 3
        for _s in sn_units_lst:
            fa_abbr_lst.extend(_s)  # Put all the FA in one list
        fa_abbr_lst = sorted(list(set(fa_abbr_lst)))

        abbr_parser = NameParserFA(lipid_class)
        elem_calc = ElemCalc()
        usr_fa_dct = {}
        for _fa_abbr in fa_abbr_lst:
            if _fa_abbr:
                _fa_info_dct = abbr_parser.get_fa_info(_fa_abbr, lipid_class, ion_mode)  # Calculate all the information for each FA
                usr_fa_dct[_fa_abbr] = _fa_info_dct

        usr_fa_df = pd.DataFrame(usr_fa_dct).T.copy()  # put all the info for the FA in a dataframe
        usr_fa_df.is_copy = False

        return usr_fa_df

    def gen_all_comb(self, lipid_class, usr_fa_df, position=False):


        fa_units_lst = self.calc_fa_df2(lipid_class, usr_fa_df)
        lipid_class_dct = NameParserFA(lipid_class).lipid_fa_dct[lipid_class]
        if lipid_class_dct[4] == 1:
            fa_comb_lst = list(itertools.product(fa_units_lst[0], ['']))
            fa_df_header_lst = lipid_class_dct[1]
            fa_df_header_lst.append('FA2')
        elif lipid_class_dct[4] == 2:
            fa_comb_lst = list(itertools.product(fa_units_lst[0], fa_units_lst[1]))
            fa_df_header_lst =lipid_class_dct[1]
            # if lipid_class == 'Cer':
            #     fa_df_header_lst = [NameParserFA().lipid_fa_dct[lipid_class][2]] + fa_df_header_lst
        elif lipid_class_dct[4] == 3:
            fa_comb_lst = list(itertools.product(fa_units_lst[0], fa_units_lst[1], fa_units_lst[2]))
            fa_df_header_lst = lipid_class_dct[1]
        elif lipid_class_dct[4] == 4:
            fa_comb_lst = list(itertools.product(fa_units_lst[0], fa_units_lst[1], fa_units_lst[2], fa_units_lst[3]))
            fa_df_header_lst = lipid_class_dct[1]
        else:
            fa_comb_lst = []
            fa_df_header_lst = []

        fa_combo_df = pd.DataFrame(data=fa_comb_lst, columns=fa_df_header_lst)

        fa_combo_df['CLASS'] = lipid_class

        # Note: a new way to combine the below section (georgia: 19.2.019)
        if 'Base' in fa_df_header_lst:
            fa_link_df = fa_combo_df
            fa_link_df.is_copy = False
            fa_link_df['DISCRETE_ABBR'] = (fa_link_df['CLASS'] + '(' + fa_link_df[fa_df_header_lst[0]] + '_' +
                                           fa_link_df[fa_df_header_lst[1]].str.strip('FA)') + ')')
            fa_link_df.sort_values(by='DISCRETE_ABBR', inplace=True)
            fa_combo_df = fa_link_df
            print ('[INFO] --> Number of predicted lipids (exact position): ', fa_combo_df.shape[0])
        else:
            fa_combo_link_df = fa_combo_df
            fa_combo_link_df.is_copy = False
            fa_combo_link_df['LINK'] = fa_combo_link_df['FA1'].str[0:2]
            fa_link_df = fa_combo_link_df[(fa_combo_link_df['LINK'] == "FA") | (fa_combo_link_df['LINK'] == "SP")]
            fa_link_df.is_copy = False
            fa_link_df.drop(['LINK', 'CLASS'], axis=1, inplace=True)
            if lipid_class_dct[2] != 0:
                fa_link_df.values.sort(kind='mergesort')
            fa_link_df['CLASS'] = lipid_class
            # Note: The header list is never len 1
            if len(fa_df_header_lst) == 1:
                fa_link_df['DISCRETE_ABBR'] = (fa_link_df['CLASS'] + '(' +
                                            fa_link_df[fa_df_header_lst[0]].str.strip('FA').str.strip('SPB') + ')')
            elif len(fa_df_header_lst) == 2:
                fa_link_df['DISCRETE_ABBR'] = (fa_link_df['CLASS'] + '(' +
                                               fa_link_df[fa_df_header_lst[0]].str.strip('FA').str.strip('SPB') + '_' +
                                               fa_link_df[fa_df_header_lst[1]].str.strip('FA') + ')')
                if lipid_class_dct[4] == 1:
                    fa_link_df['DISCRETE_ABBR'] = fa_link_df['DISCRETE_ABBR'].str.replace(r'_', '')
            elif len(fa_df_header_lst) == 3:
                fa_link_df['DISCRETE_ABBR'] = (fa_link_df['CLASS'] + '(' +
                                               fa_link_df[fa_df_header_lst[0]].str.strip('FA') + '_' +
                                               fa_link_df[fa_df_header_lst[1]].str.strip('FA') + '_' +
                                               fa_link_df[fa_df_header_lst[2]].str.strip('FA') + ')')
            elif len(fa_df_header_lst) == 4:
                fa_combo_df['DISCRETE_ABBR'] = (fa_combo_df['CLASS'] + '(' +
                                                fa_combo_df[fa_df_header_lst[0]].str.strip('FA') + '_' +
                                                fa_combo_df[fa_df_header_lst[1]].str.strip('FA') + '_' +
                                                fa_combo_df[fa_df_header_lst[2]].str.strip('FA') + '_' +
                                                fa_combo_df[fa_df_header_lst[3]].str.strip('FA') + ')')

            if lipid_class_dct[5]:
                op_link_df = fa_combo_link_df[(fa_combo_link_df['LINK'] == 'O-') | (fa_combo_link_df['LINK'] == 'P-')]
                if not op_link_df.empty:
                    op_link_df.is_copy = False
                    op_link_df.drop(['LINK'], axis = 1 , inplace=True)
                    if len(fa_df_header_lst) == 1:
                        op_link_df['DISCRETE_ABBR'] = (op_link_df['CLASS'] + '(' +
                                                       op_link_df[fa_df_header_lst[0]].str.strip('FA') + ')')
                    elif len(fa_df_header_lst) == 2:
                        op_link_df['DISCRETE_ABBR'] = (op_link_df['CLASS'] + '(' +
                                                       op_link_df[fa_df_header_lst[0]].str.strip('FA') + '_' +
                                                       op_link_df[fa_df_header_lst[1]].str.strip('FA') + ')')
                        if lipid_class_dct[4] == 1:
                            op_link_df['DISCRETE_ABBR'] = op_link_df['DISCRETE_ABBR'].str.replace(r'_', '')
                    elif len(fa_df_header_lst) == 3:
                        op_link_df['DISCRETE_ABBR'] = (op_link_df['CLASS'] + '(' +
                                                       op_link_df[fa_df_header_lst[0]].str.strip('FA') + '_' +
                                                       op_link_df[fa_df_header_lst[1]].str.strip('FA') + '_' +
                                                       op_link_df[fa_df_header_lst[2]].str.strip('FA') + ')')
                    elif len(fa_df_header_lst) == 4:
                        op_link_df['DISCRETE_ABBR'] = (op_link_df['CLASS'] + '(' +
                                                       op_link_df[fa_df_header_lst[0]].str.strip('FA') + '_' +
                                                       op_link_df[fa_df_header_lst[1]].str.strip('FA') + '_' +
                                                       op_link_df[fa_df_header_lst[2]].str.strip('FA') + '_' +
                                                       op_link_df[fa_df_header_lst[3]].str.strip('FA') + ')')
                    # TODO to the same for all the length( can also go as a function) ( georgia: 19.2.2019)
                    fa_combo_df = fa_link_df.append((op_link_df))
                else:
                    fa_combo_df = fa_link_df
            else:
                fa_combo_df = fa_link_df
            print('[INFO] --> Number of predicted lipids (exact position): ', fa_combo_df.shape[0])

        if position is False:
            print('[INFO] --> Use discrete form for identification ...')
            fa_combo_lite_df = fa_combo_df.drop_duplicates(subset=['DISCRETE_ABBR'], keep='first')
            print('[INFO] --> Number of predicted lipids (discrete form): ', fa_combo_lite_df.shape[0])
        else:
            fa_combo_lite_df = fa_combo_df

        fa_combo_lite_df.is_copy = False
        fa_combo_lite_df['idx'] = fa_combo_lite_df['DISCRETE_ABBR']
        fa_combo_lite_df.set_index('idx', drop=True, inplace=True)

        lipid_comb_dct = fa_combo_lite_df.to_dict(orient='index')

        return lipid_comb_dct


    def calc_fragments(self, lipid_dct, charge='', ms2_ppm=100):

        # m_formula = lipid_dct['FORMULA']
        abbr_name= NameParserFA(lipid_dct['CLASS'])
        # was usuful when they were the Base fragment in a separate entry in the list
        # if NameParserFA(lipid_dct['CLASS']).lipid_fa_dct[lipid_dct['CLASS']][2] is not None:
        #     fa_header_list = [NameParserFA(lipid_dct['CLASS']).lipid_fa_dct[lipid_dct['CLASS']][2]] + fa_header_list
        m_exactmass = lipid_dct['EXACTMASS']
        m_class = lipid_dct['CLASS']
        fa_regexp = re.compile(r'.*\((FA.*)\).*')
        # Note: need to check if this section is neccessary otherwise should be removed
        # Note: this section need to improve to give the correct water losses for the different ceramides base
        # Future
        lipid_frag_mode = abbr_name.lipid_mode_dct[m_class][charge]
        lipid_nl_frag = abbr_name.loss_dct

        if lipid_frag_mode[3] is not None:
            #print (NameParserFA(lipid_dct['CLASS']).lipid_mode_dct[m_class][charge][3])
            for frag in lipid_frag_mode[3]:
                fa_check = fa_regexp.match(frag)
                if fa_check:
                    fa_list_info = abbr_name.lipid_fa_dct[m_class][1]
                    for fa_i in fa_list_info:
                        nl_fr_name = '{f}_EXACTMASS'.format(f=fa_i)
                        nl_frag = lipid_dct[nl_fr_name]
                        m_exact_nl = m_exactmass - nl_frag
                        frag_n = frag.replace('X', lipid_frag_mode[4]).replace('FA', lipid_dct['{f}_ABBR'.format(f=fa_i)])
                        frag_h = frag.replace('X', lipid_frag_mode[4]).replace('FA',fa_i)

                        lipid_dct['{f}_ABBR'.format(f=frag_h)] = frag_n
                        lipid_dct['{f}_MZ'.format(f=frag_h)] = round(
                            m_exact_nl + lipid_nl_frag[frag], 6)
                        low_v = round(
                            ppm_window_para(lipid_dct['{f}_MZ'.format(f=frag_h)], ms2_ppm * -1), 6)
                        high_v = round(ppm_window_para(lipid_dct['{f}_MZ'.format(f=frag_h)], ms2_ppm), 6)
                        lipid_dct['{f}_Q'.format(f=frag_h)] = (low_v.astype(str) + ' <= mz <= '
                                                             + high_v.astype(str))
                else:
                    frag_n = frag.replace('X', lipid_frag_mode[4])
                    lipid_dct['{f}_ABBR'.format(f=frag_n)] = frag_n
                    lipid_dct['{f}_MZ'.format(f=frag_n)] = round(
                        m_exactmass + lipid_nl_frag[frag], 6)
                    low_v = round(
                            ppm_window_para(lipid_dct['{f}_MZ'.format(f=frag_n)], ms2_ppm * -1), 6)
                    high_v = round(ppm_window_para(lipid_dct['{f}_MZ'.format(f=frag_n)], ms2_ppm), 6)
                    lipid_dct['{f}_Q'.format(f=frag_n)] = (low_v.astype(str) + ' <= mz <= '
                                                                                + high_v.astype(str))
        return lipid_dct

    # Note: Get also the usr_fa_df (georgia: 20.2.2019)
    def compose_lipid(self, param_dct, usr_fa_info_df, ms2_ppm=100):

        lipid_class = param_dct['lipid_class']
        lipid_charge = param_dct['charge_mode']

        if param_dct['exact_position'] == 'TRUE':
            position_set = True
        else:
            position_set = False

        usr_fa_df = pd.read_excel(param_dct['fa_whitelist'])
        usr_fa_df = usr_fa_df.fillna(value='F')

        tmp_columns = usr_fa_df.columns.tolist()

        usr_fa_df.columns = usr_fa_df.columns.str.upper()
        # COMPINE (georgia):  the below section was rearrange to one 9.01.19

        if lipid_class in self.all_lipid_class_list:
            if lipid_class in tmp_columns:
                pass
            elif 'PL' in tmp_columns or 'LPL' in tmp_columns:
                pass
            else:
                return False
        else:
            return False

        print('[INFO] --> FA white list loaded ...')
        lipid_comb_dct = self.gen_all_comb(lipid_class, usr_fa_df, position_set)

        lipid_info_dct = {}

        elem_calc = ElemCalc()
        lipid_class_dct = NameParserFA(lipid_class).lipid_fa_dct[lipid_class]

        fa_header_lst = lipid_class_dct[1]

        for _lipid in list(lipid_comb_dct.keys()):
            _lipid_dct = lipid_comb_dct[_lipid]

            _lipid_dct['M_DB'] = 0
            _lipid_dct['M_C'] = 0
            for f_k in fa_header_lst:
                _fa_lipid = usr_fa_info_df.loc[_lipid_dct[f_k]]
                for _fa_k in list(_fa_lipid.keys()):
                    # Note: Use a filter so it would not take values that are NAN (georgia: 21.2.2019)
                    if _fa_lipid[_fa_k] is not np.nan:
                        _lipid_dct[str(f_k) + '_' + str(_fa_k)] = _fa_lipid[_fa_k]
                _lipid_dct['M_DB'] = int(_lipid_dct['M_DB']) + int(_fa_lipid['DB'])
                _lipid_dct['M_C'] = int(_lipid_dct['M_C']) + int(_fa_lipid['C'])

            # Note: thing a way to a more unify way (georgia: 20.2.2019)
            #print (usr_fa_info_df.loc[_lipid_dct[fa_header_lst[0]]])
            if usr_fa_info_df.loc[_lipid_dct[fa_header_lst[0]]]['LINK'] in ['FA', 'A']:
                lipid_bulk_str = '{pl}({c}:{db})'.format(pl=lipid_class,
                                                         c=_lipid_dct['M_C'],
                                                         db=_lipid_dct['M_DB'])
            elif usr_fa_info_df.loc[_lipid_dct[fa_header_lst[0]]]['LINK'] in ['SPB']:
                # Comment: This will not work in the case of alpha hydroxylated FA
                lipid_bulk_str = '{pl}({c}:{db};{o})'.format(pl=lipid_class,
                                                             c=_lipid_dct['M_C'],
                                                             db=_lipid_dct['M_DB'],
                                                             o= _lipid_dct['FA1_O'])
            else:
                lipid_bulk_str = '{pl}({lk}{c}:{db})'.format(pl=lipid_class,
                                                             lk= _lipid_dct[fa_header_lst[0]+'_'+'LINK'],
                                                         c=_lipid_dct['M_C'],
                                                         db=_lipid_dct['M_DB'])

            _lipid_dct['BULK_ABBR'] = lipid_bulk_str

            _lipid_formula, _lipid_elem_dct = elem_calc.get_formula(lipid_bulk_str)

            _lipid_dct['FORMULA'] = _lipid_formula
            _lipid_dct['EXACTMASS'] = elem_calc.get_exactmass(_lipid_elem_dct)
            # Note: probably this is not needed for further so could be skip the df (georgia: 20.2.2019)
            for _elem_k in list(_lipid_elem_dct.keys()):
                _lipid_dct['M_' + _elem_k] = _lipid_elem_dct[_elem_k]

            # charged
            _chg_lipid_formula, _chg_lipid_elem_dct = elem_calc.get_formula(lipid_bulk_str, charge=lipid_charge)
            _lipid_dct[lipid_charge + '_FORMULA'] = _chg_lipid_formula
            _lipid_dct[lipid_charge + '_MZ'] = elem_calc.get_exactmass(_chg_lipid_elem_dct)

            _lipid_dct = self.calc_fragments(_lipid_dct, charge=lipid_charge, ms2_ppm=ms2_ppm)

            # Question: for what reason is needed this table (georgia: 13.2.2019)
            lipid_info_dct[_lipid] = _lipid_dct
            del _lipid_dct

        # Question: how we create lipid_comb_dct (georgia: 13.2.2019)
        lipid_master_df = pd.DataFrame(lipid_comb_dct).T
        lipid_master_df.reset_index(drop=True, inplace=True)

        return lipid_master_df


if __name__ == '__main__':
    # TODO (georgia.angelidou@uni-leipzig.de): Upgrade to ceramide lipid composition
    fa_lst_file = r'../ConfigurationFiles/1-FA_Whitelist_cER_v3.xlsx'
    # fa_lst_file = r'../ConfigurationFiles/1-FA_Whitelist_TG-DG.xlsx'

    # Note:
    # exact position means to consider the position from the FA white list that the user give but,
    # in the case that the user define 2 different FA for both positions then:
    # When it is false it will give only one option
    # and when it is TRUE to give both combinations that these 2 FA an make (in case of phospholipids)

    # usr_param_dct = {'fa_whitelist': fa_lst_file, 'lipid_class': 'Cer', 'charge_mode': '[M+H]+',
    #                  'exact_position': 'FALSE'}

    usr_param_dct = {'fa_whitelist': fa_lst_file, 'lipid_class': 'Cer', 'charge_mode': '[M+H]+',
                     'exact_position': 'FALSE'}
    # usr_param_dct = {'fa_whitelist': fa_lst_file, 'lipid_class': 'TG', 'charge_mode': '[M+NH4]+',
    #                  'exact_position': 'FALSE'}

    composer = LipidComposer()
    calc_fa_df = composer.calc_fa_query(usr_param_dct['lipid_class'], usr_param_dct['charge_mode'],
                                        fa_lst_file, ms2_ppm=50)
    usr_lipid_master_df = composer.compose_lipid(usr_param_dct,  calc_fa_df,  ms2_ppm=30)
    print('[INFO] --> Lipid Master Table generated...')
    master_csv = r'../Temp/LipidMaster_Whitelist_%s.csv' % usr_param_dct['lipid_class']
    usr_lipid_master_df.to_csv(master_csv)
    # master_csv = r'../Temp/LipidMaster_Whitelist_TG_ML.csv'

    fa_csv = r'../Temp/LipidMaster_FAlist2_%s.csv' % usr_param_dct['lipid_class']

    # calc_fa_df = composer.calc_fa_query(usr_param_dct['lipid_class'], usr_param_dct['charge_mode'],
    #                                     fa_lst_file, ms2_ppm=50)

    if calc_fa_df is False:
        print('[ERROR] !!! Failed to generate FA info table ...\n')

    print(calc_fa_df.head())

    calc_fa_df.to_csv(fa_csv)
    print('[INFO] --> Finished...')
