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

import re

try:
    from LibLipidHunter.ParallelFunc import ppm_window_para
    from LibLipidHunter.AbbrElemCalc import ElemCalc
except ImportError:
    from ParallelFunc import ppm_window_para
    from AbbrElemCalc import ElemCalc


class NameParserFA:
    # TODO (georgia): need to make this more general
    def __init__(self, lclass = 'PC'):
        self.elem_dct = {'H': [1.0078250321, 0.999885],
                         'D': [2.0141017780, 0.0001157],
                         'C': [12.0, 0.9893],
                         'N': [14.0030740052, 0.99632],
                         'O': [15.9949146221, 0.99757],
                         'Na': [22.98976967, 1.0],
                         'P': [30.97376151, 1.0],
                         'S': [31.97207069, 0.9493],
                         'K': [38.9637069, 0.932581]}

        # The regular expression discrebe all the possible name abbreviation for the free lipids including ceramide base
        # First Section: ([a-zA-z\-]{1,2}) -> for the following options ('FA', 'O-', 'P-', 'd', 'm', 't')
        # The abbreviation ('FA', 'O-', 'P-') indicate if there normal fatty acid or plasmalogen
        # The abbreviation ('d', 'm', 't') indicate how many hydroxy groups the ceramide base(sphingoninnes) have
        # Second Section: {1,2})(\d{1,2})(:)(\d) -> For the total number of C and DB of the lipid chains
        self.fa_rgx = re.compile(r'([a-zA-z\-]{1,3})(\d{1,2})(:)(\d)(;)*(\d)*')

        neg_gen = ['[FA-H]-', '[FA-H2O]', '[FA-H2O-H]-']
        neg_gen2 = ['[x-H]-', '[x-H2O-H]-']
        pc_neg = ['[x-CH3]-', '[x-H2O-CH3]-']
        # Note: MG fragment can be move to the pos_gen or to pos_gen2 (georgia: 22.2.2019)
        pos_gen = ['[FA-H2O+H]+', '[FA-H2O]']
        pos_gen2 = ['[x-H2O+H]+']
        pos_sod_gen_the2 = ['[X-(FA)+Na]+', '[X-(FA-H+Na)+H]+']
        pos_gen_the2 = ['[X-(FA)+H]+', '[X-(FA-H2O)+H]+']
        sodium_spec = ['[FA-H+Na]']
        cer_base_spec = ['[LCB+H]+','[LCB-H2O+H]+', '[LCB-2xH2O+H]+', '[LCB-H2O-CH2O+H]+', '[LCB-3xH2O+H]+', '[LCB-2xH2O-CH2O+H]+']
        cer_base_spec2 = ['[X-H2O+H]+', '[X-2xH2O+H]+', '[X-H2O-CH2O+H]+', '[X-3xH2O+H]+', '[X-2xH2O-CH2O+H]+']
        cer_fa_spec = ['[FA-OH+NH3]+', '[FA-OH+C2H3N]+']

        self.loss_dct = {'[FA-H]-': -self.elem_dct['H'][0],
                         '[FA-H2O]': -(self.elem_dct['H'][0]*2 + self.elem_dct['O'][0]),
                         '[FA-H2O-H]-': -(self.elem_dct['H'][0]*3 + self.elem_dct['O'][0]),
                         '[FA-H2O+H]+': -(self.elem_dct['H'][0] + self.elem_dct['O'][0]),
                         '[FA-H+Na]': self.elem_dct['Na'][0] - self.elem_dct['H'][0],
                         '[LCB+H]+': self.elem_dct['H'][0],
                         '[LCB-H2O+H]+': -(self.elem_dct['H'][0] + self.elem_dct['O'][0]),
                         '[LCB-2xH2O+H]+': -(self.elem_dct['H'][0]*3 + self.elem_dct['O'][0]*2),
                         '[LCB-H2O-CH2O+H]+': -(self.elem_dct['H'][0]*3 + self.elem_dct['O'][0]*2 +
                                                self.elem_dct['C'][0]),
                         '[LCB-3xH2O+H]+': -(self.elem_dct['H'][0] * 5 + self.elem_dct['O'][0] * 3),
                         '[LCB-2xH2O-CH2O+H]+': -(
                                     self.elem_dct['H'][0] * 5 + self.elem_dct['O'][0] * 3 + self.elem_dct['C'][0]),
                         '[X+H]+': self.elem_dct['H'][0],
                         '[X-H2O+H]+': -(self.elem_dct['H'][0] + self.elem_dct['O'][0]),
                         '[X-2xH2O+H]+': -(self.elem_dct['H'][0] * 3 + self.elem_dct['O'][0] * 2),
                         '[X-H2O-CH2O+H]+': -(self.elem_dct['H'][0] * 3 + self.elem_dct['O'][0] * 2 +
                                                self.elem_dct['C'][0]),
                         '[X-3xH2O+H]+': -(self.elem_dct['H'][0] * 5 + self.elem_dct['O'][0] * 3),
                         '[X-2xH2O-CH2O+H]+': -(
                                 self.elem_dct['H'][0] * 5 + self.elem_dct['O'][0] * 3 + self.elem_dct['C'][0]),
                         '[FA-OH+NH3]+': self.elem_dct['N'][0] + self.elem_dct['H'][0]*2 - self.elem_dct['O'][0],
                         '[FA-OH+C2H3N]+': (self.elem_dct['C'][0] + self.elem_dct['N'][0] + self.elem_dct['H'][0]*2 -
                                            self.elem_dct['O'][0]),
                         '[x-H]-': self.calc_fa_mass(ElemCalc().lipid_hg_elem_dct[lclass])
                                   + self.calc_fa_mass(ElemCalc().glycerol_bone_elem_dct)
                                   + self.elem_dct['H'][0] + self.elem_dct['O'][0],
                         '[x-H2O-H]-': self.calc_fa_mass(ElemCalc().lipid_hg_elem_dct[lclass])
                                   + self.calc_fa_mass(ElemCalc().glycerol_bone_elem_dct)
                                   - self.elem_dct['H'][0],
                         '[x-CH3]-': self.calc_fa_mass(ElemCalc().lipid_hg_elem_dct[lclass])
                                   + self.calc_fa_mass(ElemCalc().glycerol_bone_elem_dct)
                                   - self.elem_dct['H'][0] + self.elem_dct['O'][0] - self.elem_dct['C'][0],
                         '[x-H2O-CH3]-': self.calc_fa_mass(ElemCalc().lipid_hg_elem_dct[lclass])
                                   + self.calc_fa_mass(ElemCalc().glycerol_bone_elem_dct) - self.elem_dct['C'][0]
                                         - self.elem_dct['H'][0]*3,
                         '[x-H2O+H]+': self.calc_fa_mass(ElemCalc().glycerol_bone_elem_dct) + self.elem_dct['H'][0]*3
                                       + self.elem_dct['O'][0],
                         '[X-(FA)+H]+': self.elem_dct['H'][0],  # will be calculate from the exact neutral mass of the lipid thats the reason we do not need to add the neutral loss of NH3 in the case of TG
                         '[X-(FA-H2O)+H]+': self.elem_dct['H'][0] * 3 + self.elem_dct['O'][0],  # Will be calculated by the exact neutral mass of the lipid molecule minus the exact mass of the fa.
                                                                                                # In this case where the fa is loss without the water we need to add it in this step to have the correct calculation
                         '[X-(FA)+Na]+': self.elem_dct['Na'][0],
                         '[X-(FA-H+Na)+H]+': self.elem_dct['H'][0]}



        self.lipid_mode_dct = {'PA': {'[M-H]-': [neg_gen, neg_gen2, 'LPA', None, None]},
                               'PC': {'[M+HCOO]-': [neg_gen, pc_neg, 'LPC', None, None],
                                      '[M+CH3COO]-': [neg_gen, pc_neg, 'LPC', None, None]},
                               'PE': {'[M-H]-': [neg_gen, neg_gen2, 'LPE', None, None]},
                               'PG': {'[M-H]-': [neg_gen, neg_gen2, 'LPG', None, None]},
                               'PI': {'[M-H]-': [neg_gen, neg_gen2, 'LPI', None, None]},
                               'PS': {'[M-H]-': [neg_gen, neg_gen2, 'LPS', None, None]},
                               'TG': {'[M+NH4]+': [pos_gen, pos_gen2, 'MG', pos_gen_the2, 'M'],
                                      '[M+H]+': [pos_gen, pos_gen2, 'MG', pos_gen_the2, 'M'],
                                      '[M+Na]+': [(pos_gen + sodium_spec), pos_gen2, 'MG', pos_sod_gen_the2 , 'M']},
                               'DG': {'[M+NH4]+': [pos_gen, pos_gen2, 'MG', None, None]},
                               'LPA': {'[M-H]-': [neg_gen, neg_gen2, 'LPA', None, None]},
                               'LPC': {'[M+HCOO]-': [neg_gen, neg_gen2, 'LPC', None, None],
                                       '[M+CH3COO]-': [neg_gen, neg_gen2, 'LPC', None, None]},
                               'LPE': {'[M-H]-': [neg_gen, neg_gen2, 'LPE', None, None]},
                               'LPG': {'[M-H]-': [neg_gen, neg_gen2, 'LPG', None, None]},
                               'LPI': {'[M-H]-': [neg_gen, neg_gen2, 'LPI', None, None]},
                               'LPS': {'[M-H]-': [neg_gen, neg_gen2, 'LPS', None, None]},
                               'Cer': {'[M+H]+': [cer_fa_spec, cer_base_spec, 'LCB', cer_base_spec2, 'M']}}
        # Note: add also the length of the list that it should have (georgia: 18.2.2019)
        # Note: in the below section we can add also the info about the link. (georgia: 18.2.2019)
        # Note: The true or false statement (last entry of the list indicate if we have O- or P- FA
        # Note: Entry at the position 3 in the list will be used as a control if should short the FA or not
        # Note (continue): 0 - indicates to do not mergesort them
        # Note: The position 4 in the define list below should indicate if the first FA1 should be sphingoid base or not.
        # Note (c): this way the program will identify if the general class is sphingolipids or not and avoid mistakes
        self.lipid_fa_dct = {'PA': ['PL', ['FA1', 'FA2'], None, None, 2, False],
                             'PC': ['PL', ['FA1', 'FA2'], None, None, 2, True],
                             'PE': ['PL', ['FA1', 'FA2'], None, None, 2, True],
                             'PG': ['PL', ['FA1', 'FA2'], None, None, 2, False],
                             'PI': ['PL', ['FA1', 'FA2'], None, None, 2, False],
                             'PS': ['PL', ['FA1', 'FA2'], None, None, 2, False],
                             'TG': ['TG', ['FA1', 'FA2', 'FA3'], None, None, 3, True],
                             'DG': ['DG', ['FA1', 'FA2'], None, None, 2, True],
                             'LPA': ['LPL', ['FA1'], 0, None, 1, False],
                             'LPC': ['LPL', ['FA1'], 0, None, 1, True],
                             'LPE': ['LPL', ['FA1'], 0, None, 1, True],
                             'LPG': ['LPL', ['FA1'], 0, None, 1, False],
                             'LPI': ['LPL', ['FA1'], 0, None, 1, False],
                             'LPS': ['LPL', ['FA1'], 0, None, 1, False],
                             'CL': ['CL', ['FA1', 'FA2', 'FA3', 'FA4'], None, None, 4, False],
                             'Cer': ['CER', ['FA1', 'FA2'], 0, 'SPB', 2, False]} # alternative need to decide which of the two wil be kept
                             # 'Cer': ['CER', ['FA1'], 'Base', ['M', 'D', 'T'], 2, False]}

        self.lipid_link_dct = {'PA': False}
        # Now create manually but consideraration to create this more dynamically
        lpl_a = 'L{l}'.format(l=lclass)
        self.lipid_fragment_info = {'[LPL-H]-': ['FA1_[LPC-CH3]-', 'FA2_[LPC-CH3]-',
                                                     'FA1_[{l}-H]-'.format(l=lpl_a), 'FA2_[{l}-H]-'.format(l=lpl_a)],
                                    '[LPL-H2O-H]-': ['FA1_[LPC-H2O-CH3]-', 'FA2_[LPC-H2O-CH3]-',
                                                     'FA1_[{l}-H2O-H]-'.format(l=lpl_a),
                                                     'FA2_[{l}-H2O-H]-'.format(l=lpl_a)],
                                    '[FA-H]-': ['FA1_[FA-H]-', 'FA2_[FA-H]-'],
                                    '[FA-H2O+H]+': ['FA1_[FA-H2O+H]+', 'FA2_[FA-H2O+H]+', 'FA3_[FA-H2O+H]+'],
                                    '[FA-H2O-H]+': ['FA1_[FA-H2O-H]-', 'FA2_[FA-H2O-H]-', 'FA3_[FA-H2O-H]-'],
                                    '[MG-H2O+H]+': ['FA1_[MG-H2O+H]+', 'FA2_[MG-H2O+H]+', 'FA3_[MG-H2O+H]+'],
                                    '[M-FA+H]+': ['[M-(FA1)+H]+', '[M-(FA2)+H]+', '[M-(FA3)+H]+'],
                                    '[M-(FA-H2O)+H]+': ['[M-(FA1-H2O)+H]+', '[M-(FA2-H2O)+H]+', '[M-(FA3-H2O)+H]+'],
                                    '[M-FA+Na]+': ['[M-(FA1)+Na]+', '[M-(FA2)+Na]+', '[M-(FA3)+Na]+'],
                                    '[M-(FA-H+Na)+H]+': ['[M-(FA1-H+Na)+H]+', '[M-(FA2-H+Na)+H]+',
                                                         '[M-(FA3-H+Na)+H]+'],
                                    '[LCB...+H]+': ['FA1_[LCB+H]+', 'FA1_[LCB-H2O+H]+', 'FA1_[LCB-2xH2O+H]+',
                                                    'FA1_[LCB-H2O-CH2O+H]+', 'FA1_[LCB-3xH2O+H]+',
                                                    'FA1_[LCB-2xH2O-CH2O+H]+'],
                                    '[M...+H]+': ['[M-H2O+H]+', '[M-2xH2O+H]+', '[M-3xH2O+H]+', '[M-H2O-CH2O+H]+',
                                                  '[M-2xH2O-CH2O+H]+'],
                                    '[FA-OH+NH3]+': ['FA2_[FA-OH+NH3]+']}
        self.lipid_output_format = {'PC': {'[M+HCOO]-':{'OBS_FA': [['[FA-H]-', 'upper left'], ['[FA-H2O-H]-', 'upper right']],
                                                        'OBS_LYSO': [['[LPL-H]-', 'upper left'], ['[LPL-H2O-H]-', 'upper right']] },
                                           '[M+CH3COO]-':{'OBS_FA': [['[FA-H]-', 'upper left']],
                                                        'OBS_LYSO': [['[LPL-H]-', 'upper left'], ['[LPL-H2O-H]-', 'upper right']]}},
                                    'PE': {'[M-H]-':{'OBS_FA': [['[FA-H]-', 'upper left']],
                                                     'OBS_LYSO': [['[LPL-H]-', 'upper left'], ['[LPL-H2O-H]-', 'upper right']]}},
                                    'PS': {'[M-H]-':{'OBS_FA': [['[FA-H]-', 'upper left']],
                                                     'OBS_LYSO': [['[LPL-H]-', 'upper left'], ['[LPL-H2O-H]-', 'upper right']]}},
                                    'PG': {'[M-H]-':{'OBS_FA': [['[FA-H]-', 'upper left']],
                                                     'OBS_LYSO': [['[LPL-H]-', 'upper left'], ['[LPL-H2O-H]-', 'upper right']]}},
                                    'PI': {'[M-H]-': {'OBS_FA': [['[FA-H]-', 'upper left']],
                                                      'OBS_LYSO': [['[LPL-H]-', '[LPI-H2O-H]-']]}},
                                    'PA': {'[M-H]-': {'OBS_FA': [['[FA-H]-', 'upper left']],
                                                      'OBS_LYSO': [['[LPL-H]-', 'upper left'], ['[LPL-H2O-H]-', 'upper right']]}},
                                    'DG': {'[M+NH4]+': {'OBS_FA': [['[FA-H2O+H]+', 'upper left']],
                                                      'OBS_LYSO': [['[MG-H2O+H]+', 'upper left']]}},
                                    'TG': {'[M+NH4]+': {'OBS_FA': [['[FA-H2O+H]+', 'upper left'], ['[MG-H2O+H]+', 'upper right']],
                                                        'OBS_LYSO': [['[M-FA+H]+', 'upper left']]},
                                           '[M+Na]+': {'OBS_FA': [['[FA-H2O+H]+', 'upper left'], ['[MG-H2O+H]+', 'upper right']],
                                                       'OBS_LYSO': [['[M-FA+Na]+', 'upper left'], ['[M-(FA-H+Na)+H]+', 'upper right']]}},
                                    'Cer': {'[M+H]+': {'OBS_FA': [['[LCB...+H]+', 'upper left'], ['[FA-OH+NH3]+', 'upper right']],
                                                       'OBS_LYSO': [['[M...+H]+', 'upper left']]}},
                                    'LPC': {'[M+HCOO]-':{'OBS_FA': [['[FA-H]-', 'upper left']],
                                                        'OBS_LYSO': []},
                                           '[M+CH3COO]-':{'OBS_FA': [['[FA-H]-', 'upper left']],
                                                        'OBS_LYSO': []}},
                                    'LPE': {'[M-H]-':{'OBS_FA': [['[FA-H]-', 'upper left']],
                                                     'OBS_LYSO': []}},
                                    'LPS': {'[M-H]-':{'OBS_FA': [['[FA-H]-', 'upper left']],
                                                     'OBS_LYSO': []}},
                                    'LPG': {'[M-H]-':{'OBS_FA': [['[FA-H]-', 'upper left']],
                                                     'OBS_LYSO': []}},
                                    'LPI': {'[M-H]-': {'OBS_FA': [['[FA-H]-', 'upper left']],
                                                      'OBS_LYSO': []}},
                                    'LPA': {'[M-H]-': {'OBS_FA': [['[FA-H]-', 'upper left']],
                                                      'OBS_LYSO': []}}}

    # Note: Current version use only once but maybe will be usefull in the future
    def calc_fa_mass(self, fa_info_dct):
        exactmass = 0
        # Note: this version causing issues with Cer since the also have an N for this reason re-arrange
        # exactmass = (fa_info_dct['C'] * self.elem_dct['C'][0] + fa_info_dct['H'] * self.elem_dct['H'][0] +
        #              fa_info_dct['O'] * self.elem_dct['O'][0])
        # Re-arrange: (georgia-23.1.2019)
        for k in fa_info_dct.keys():
            if k in self.elem_dct.keys():
                exactmass = exactmass + (fa_info_dct[k] * self.elem_dct[k][0])
        return round(exactmass, 6)

    def get_fa_formula(self, fa_str):

        fa_info_dct = {}
        # fa_info_lst = self.decode_fa(fa_str)
        fa_info_lst = re.match(self.fa_rgx, fa_str).groups()
        if fa_info_lst[0] == 'FA':
            fa_info_dct['ABBR'] = fa_str
            fa_info_dct['LINK'] = 'FA'
            # Note: Below section can be skip and just call get_neutral_elem
            fa_info_dct['C'] = int(fa_info_lst[1])
            fa_info_dct['H'] = 2 * int(fa_info_lst[1]) - 2 * int(fa_info_lst[3])
            fa_info_dct['O'] = 2
            fa_info_dct['DB'] = int(fa_info_lst[3])
            fa_info_dct['FORMULA'] = 'C{num_c}H{num_h}O{num_o}'.format(num_c=fa_info_dct['C'], num_h=fa_info_dct['H'],
                                                                     num_o=fa_info_dct['O'])
        elif fa_info_lst[0] == 'O-':
            fa_info_dct['ABBR'] = fa_str
            fa_info_dct['LINK'] = 'O-'
            fa_info_dct['C'] = int(fa_info_lst[1])
            fa_info_dct['H'] = 2 * int(fa_info_lst[1]) + 2 - 2 * int(fa_info_lst[3])
            fa_info_dct['O'] = 1
            fa_info_dct['DB'] = int(fa_info_lst[3])
            fa_info_dct['FORMULA'] = 'C{num_c}H{num_h}O'.format(num_c=fa_info_dct['C'], num_h=fa_info_dct['H'])
        elif fa_info_lst[0] == 'P-':
            fa_info_dct['ABBR'] = fa_str
            fa_info_dct['LINK'] = 'P-'
            fa_info_dct['C'] = int(fa_info_lst[1])
            fa_info_dct['H'] = 2 * int(fa_info_lst[1]) - 2 * int(fa_info_lst[3])
            fa_info_dct['O'] = 1
            fa_info_dct['DB'] = int(fa_info_lst[3])
            fa_info_dct['FORMULA'] = 'C{num_c}H{num_h}O'.format(num_c=fa_info_dct['C'], num_h=fa_info_dct['H'])
        elif fa_info_lst[0] == 'SPB':
            fa_info_dct['LINK'] = 'SPB'
            fa_info_dct['ABBR'] = fa_str
            fa_info_dct['C'] = int(fa_info_lst[1])
            fa_info_dct['H'] = 2 * int(fa_info_lst[1]) - 2 * int(fa_info_lst[3]) + 3
            fa_info_dct['O'] = int(fa_info_lst[5])
            fa_info_dct['N'] = 1
            fa_info_dct['DB'] = fa_info_lst[3]
            fa_info_dct['FORMULA'] = 'C{num_c}H{num_h}O{num_o}N'.format(num_c=fa_info_dct['C'], num_h=fa_info_dct['H'],
                                                                       num_o=fa_info_dct['O'])
        return fa_info_dct

    # Note: More general way to define the fragments. Also define specific and not. (georgia: 23.1.2019)
    #Change: Values for HIGH and LOW can be removed. Check if are need or not (georgia: 15.2.2019)
    # TODO: check if they cause an issue or not. (low_v & high_v) (georgia: 15.2.2019)
    def frag_definition(self, lclass, mode, fa_info , ms2_ppm = 20):
        #print ('lets see what we are going to do here')
        exactmass = self.calc_fa_mass(fa_info)
        fa_info['EXACTMASS'] = exactmass
        low_v = ''
        high_v = ''

        if fa_info['LINK'] == 'SPB':
            #print ('welll')
            fragment_list = self.lipid_mode_dct[lclass][mode][1]
            if fa_info['O'] == 1:
                sub_fragment_list = fragment_list[0:2]
            elif fa_info['O'] == 2:
                sub_fragment_list = fragment_list[0:4]
            elif fa_info['O'] == 3:
                sub_fragment_list = fragment_list[0:3] + fragment_list[4:6]
            c_name = self.lipid_mode_dct[lclass][mode][2]
            for frag in sub_fragment_list:
                fa_info['{f}_MZ'.format(f=frag)]= round(exactmass + self.loss_dct[frag], 6)
                fa_info['{f}_ABBR'.format(f=frag)] = frag.replace('LCB', '{c}({l})'.format(c=c_name, l=fa_info['ABBR'].replace('SPB', '')))
                low_v = round(ppm_window_para(fa_info['{f}_MZ'.format(f=frag)], ms2_ppm * -1), 6)
                high_v = round(ppm_window_para(fa_info['{f}_MZ'.format(f=frag)], ms2_ppm), 6)
                fa_info['{f}_Q'.format(f=frag)] = (low_v.astype(str) + ' <= mz <= ' + high_v.astype(str))


        elif fa_info['LINK'] in ['FA', 'O-', 'P-']:
                fragment_list = self.lipid_mode_dct[lclass][mode][0]
                for frag in fragment_list:
                    fa_info['{f}_MZ'.format(f=frag)] = round(exactmass + self.loss_dct[frag], 6)
                    fa_info['{f}_ABBR'.format(f=frag)] = frag.replace('FA', fa_info['ABBR'])
                    low_v = round(ppm_window_para(fa_info['{f}_MZ'.format(f=frag)], ms2_ppm * -1), 6)
                    high_v = round(ppm_window_para(fa_info['{f}_MZ'.format(f=frag)], ms2_ppm), 6)
                    fa_info['{f}_Q'.format(f=frag)] = (low_v.astype(str) + ' <= mz <= ' + high_v.astype(str))

                if lclass != 'Cer':
                    fragment_list2 = self.lipid_mode_dct[lclass][mode][1]
                    fragm2_name = "{l}({f})".format(l=self.lipid_mode_dct[lclass][mode][2], f=fa_info['ABBR'].replace('FA', ''))

                    for frag in fragment_list2:
                        fragm_name = frag.replace('x', self.lipid_mode_dct[lclass][mode][2])
                        fa_info['{f}_MZ'.format(f=fragm_name)] = round(exactmass + self.loss_dct[frag], 6)
                        fa_info['{f}_ABBR'.format(f=fragm_name)] = frag.replace('x', fragm2_name)
                        low_v = round(ppm_window_para(fa_info['{f}_MZ'.format(f=fragm_name)],
                                                                               ms2_ppm * -1), 6)
                        high_v = round(ppm_window_para(fa_info['{f}_MZ'.format(f=fragm_name)], ms2_ppm), 6)
                        fa_info['{f}_Q'.format(f=fragm_name)] = (low_v.astype(str)) + ' <= mz <= ' + (high_v.astype(str))

        return (fa_info)
    # Note: my function so i will be able to test a theory about who to solve a problem
    def get_fa_info(self, fa_str, lipid_class, ion_mode, ms2_ppm = 20):

        #print (fa_str)
        fa_info_dct = self.get_fa_formula(fa_str)
        fa_info_dct = self.frag_definition(lipid_class, ion_mode, fa_info_dct)


        return fa_info_dct


if __name__ == '__main__':


    abbr_decoder = NameParserFA()
    lipidC = 'Cer'
    ionM = '[M+H]+'
    abbr_lst = ['SPB18:1;1', 'SPB18:1;2',  'SPB18:1;3', 'FA16:0', 'FA18:0']
    for abbr in abbr_lst:
        x = abbr_decoder.get_fa_info(abbr, lipidC, ionM)
        print(x)
    # Note: Alternative way to test the function
    # from test.test_LipidNomenclature import test_get_fa_info
    #
    # test_get_fa_info()

