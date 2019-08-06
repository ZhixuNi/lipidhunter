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

        # Compine: (georgia_10.1.2019): Below section compine in one line
        # Delet: (georgia_10.1.2019)
        # self.fa_rgx = re.compile(r'(\w{1,2})(\d{1,2})(:)(\d)')
        # self.o_rgx = re.compile(r'(O-)(\d{1,2})(:)(\d)')
        # self.p_rgx = re.compile(r'(P-)(\d{1,2})(:)(\d)')
        # self.fa_rgx_lst = [self.fa_rgx, self.o_rgx, self.p_rgx]

        # The regular expression discrebe all the possible name abbreviation for the free lipids including ceramide base
        # First Section: ([a-zA-z\-]{1,2}) -> for the following options ('FA', 'O-', 'P-', 'd', 'm', 't')
        # The abbreviation ('FA', 'O-', 'P-') indicate if there normal fatty acid or plasmalogen
        # The abbreviation ('d', 'm', 't') indicate how many hydroxy groups the ceramide base(sphingoninnes) have
        # Second Section: {1,2})(\d{1,2})(:)(\d) -> For the total number of C and DB of the lipid chains
        self.fa_rgx = re.compile(r'([a-zA-z\-]{1,2})(\d{1,2})(:)(\d)')

        neg_gen = ['[FA-H]-', '[FA-H2O]', '[FA-H2O-H]-']
        neg_gen2 = ['[x-H]-', '[x-H2O-H]-']
        pc_neg = ['[x-CH3]-', '[x-H2O-CH3]-']
        # Note: MG fragment can be move to the pos_gen or to pos_gen2 (georgia: 22.2.2019)
        pos_gen = ['[FA-H2O+H]+', '[FA-H2O]']
        # pos_gen = ['[FA-H2O+H]+', '[FA-H2O]', '[MG(FA)-H2O+H]+']
        pos_gen2 = ['[x-H2O+H]+']
        pos_sod_gen_the = ['[M-(FA)+Na]+', '[M-(FA-H+Na)+H]+']
        pos_gen_the = ['[M-(FA)+H]+', '[M-(FA-H2O)+H]+']
        sodium_spec = ['[FA-H+Na]']
        sodium_spec2 = ['[x-H+Na]']
        cer_base_spec = ['[LCB-H2O+H]+', '[LCB-2xH2O+H]+', '[LCB-H2O-CH2O+H]+']
        cer_base_spec2 = ['[x-H2O+H]+', '[x-2xH2O+H]+', '[x-H2O-CH2O+H]+']
        cer_fa_spec = ['[FA-OH+NH3]+', '[FA-OH+C2H3N]+']
        cer_fa_spec2 = ['[x-OH+NH3]+', '[x-OH+C2H3N]+']


        self.loss_dct = {'[FA-H]-': -self.elem_dct['H'][0],
                         '[FA-H2O]': -(self.elem_dct['H'][0]*2 + self.elem_dct['O'][0]),
                         '[FA-H2O-H]-': -(self.elem_dct['H'][0]*3 + self.elem_dct['O'][0]),
                         '[FA-H2O+H]+': -(self.elem_dct['H'][0] + self.elem_dct['O'][0]),
                         '[FA-H+Na]': self.elem_dct['Na'][0] - self.elem_dct['H'][0],
                         '[LCB-H2O+H]+': -(self.elem_dct['H'][0] + self.elem_dct['O'][0]),
                         '[LCB-2xH2O+H]+': -(self.elem_dct['H'][0]*3 + self.elem_dct['O'][0]*2),
                         '[LCB-H2O-CH2O+H]+': -(self.elem_dct['H'][0]*3 + self.elem_dct['O'][0]*2 + self.elem_dct['C'][0]),
                         '[FA-OH+NH3]+': self.elem_dct['N'][0] + self.elem_dct['H'][0]*2 - self.elem_dct['O'][0],
                         '[FA-OH+C2H3N]+': (self.elem_dct['C'][0] + self.elem_dct['N'][0] + self.elem_dct['H'][0]*2 -
                                            self.elem_dct['O'][0]),
                         '[MG(FA)-H2O+H]+': self.calc_fa_mass(ElemCalc().lipid_hg_elem_dct[lclass])
                                   + self.calc_fa_mass(ElemCalc().glycerol_bone_elem_dct)
                                   - self.elem_dct['H'][0],
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
                         '[M-(FA)+H]+': self.elem_dct['H'][0],
                         '[M-(FA-H2O)+H]+': self.elem_dct['H'][0]*3 + self.elem_dct['O'][0]}
        # Note: the reason that i cannot use the above is because of the adiotion of glycerol bone or head (georgia: 21.2.2019)
        # Note: Need to find a more unify way who to do this (georgia: 21.2.2019)
        # Note: one solution is to put the above in a list where the first cell is the loss and the other the fragments info (georgia: 21.2.2019)

        self.loss_dct2 = {'[FA-H]-': [-self.elem_dct['H'][0], ''],
                          '[FA-H2O]': [-(self.elem_dct['H'][0] * 2 + self.elem_dct['O'][0]), ''],
                          '[FA-H2O-H]-': [-(self.elem_dct['H'][0] * 3 + self.elem_dct['O'][0]), ''],
                          '[FA-H2O+H]+': [-(self.elem_dct['H'][0] + self.elem_dct['O'][0]), ''],
                          '[FA-H+Na]': [self.elem_dct['Na'][0] - self.elem_dct['H'][0], ''],
                          '[LCB-H2O+H]+': [-(self.elem_dct['H'][0] + self.elem_dct['O'][0]), ''],
                          '[LCB-2xH2O+H]+': [-(self.elem_dct['H'][0] * 3 + self.elem_dct['O'][0] * 2), ''],
                          '[LCB-H2O-CH2O+H]+': [-(
                                     self.elem_dct['H'][0] * 3 + self.elem_dct['O'][0] * 2 + self.elem_dct['C'][0]), ''],
                          '[FA-OH+NH3]+': [self.elem_dct['N'][0] + self.elem_dct['H'][0] * 2 - self.elem_dct['O'][0], ''],
                          '[FA-OH+C2H3N]+': [(self.elem_dct['C'][0] + self.elem_dct['N'][0] + self.elem_dct['H'][0] * 2 -
                                            self.elem_dct['O'][0]), ''],
                          '[MG(FA)-H2O+H]+': [ -self.elem_dct['H'][0], self.calc_fa_mass(ElemCalc().lipid_hg_elem_dct[lclass])
                                             + self.calc_fa_mass(ElemCalc().glycerol_bone_elem_dct)],
                          '[x-H]-': [self.elem_dct['H'][0] + self.elem_dct['O'][0], self.calc_fa_mass(ElemCalc().lipid_hg_elem_dct[lclass])
                                   + self.calc_fa_mass(ElemCalc().glycerol_bone_elem_dct)],
                          '[x-H2O-H]-': [ -self.elem_dct['H'][0], self.calc_fa_mass(ElemCalc().lipid_hg_elem_dct[lclass])
                                       + self.calc_fa_mass(ElemCalc().glycerol_bone_elem_dct)],
                          '[x-CH3]-': [ -self.elem_dct['H'][0] + self.elem_dct['O'][0] - self.elem_dct['C'][0],
                                       self.calc_fa_mass(ElemCalc().lipid_hg_elem_dct[lclass])
                                     + self.calc_fa_mass(ElemCalc().glycerol_bone_elem_dct)],
                          '[x-H2O-CH3]-': [ -self.elem_dct['C'][0] - self.elem_dct['H'][0] * 3,
                                           self.calc_fa_mass(ElemCalc().lipid_hg_elem_dct[lclass])
                                         + self.calc_fa_mass(ElemCalc().glycerol_bone_elem_dct)],
                          '[x-H2O+H]+': [self.elem_dct['H'][0] * 3 + self.elem_dct['O'][0],
                                        self.calc_fa_mass(ElemCalc().glycerol_bone_elem_dct)]}




        self.lipid_mode_dct = {'PA': {'[M-H]-': [neg_gen, neg_gen2, 'LPA', neg_gen2, 'LPL']},
                               'PC': {'[M+HCOO]-': [neg_gen, pc_neg, 'LPC', neg_gen2, 'LPL'],
                                      '[M+CH3COO]-': [neg_gen, pc_neg, 'LPC', neg_gen2, 'LPL']},
                               'PE': {'[M-H]-': [neg_gen, neg_gen2, 'LPE', neg_gen2, 'LPL']},
                               'PG': {'[M-H]-': [neg_gen, neg_gen2, 'LPG', neg_gen2, 'LPL']},
                               'PI': {'[M-H]-': [neg_gen, neg_gen2, 'LPI', neg_gen2, 'LPL']},
                               'PS': {'[M-H]-': [neg_gen, neg_gen2, 'LPS', neg_gen2, 'LPL']},
                               'TG': {'[M+NH4]+': [pos_gen, pos_gen2, 'MG', pos_gen_the, None],
                                      '[M+H]+': [pos_gen, pos_gen2, 'MG', pos_gen_the, None],
                                      '[M+Na]+': [(pos_gen + sodium_spec), pos_gen2, 'MG', pos_sod_gen_the , None]},
                               'DG': {'[M+NH4]+': [pos_gen, pos_gen2, 'MG', None, None]},
                               'LPA': {'[M-H]-': [neg_gen, neg_gen2, lclass, None, None]},
                               'LPC': {'[M+HCOO]-': [neg_gen, neg_gen2, lclass, None, None],
                                       '[M+CH3COO]-': [neg_gen, neg_gen2, lclass, None, None]},
                               'LPE': {'[M-H]-': [neg_gen, neg_gen2, lclass, None, None]},
                               'LPG': {'[M-H]-': [neg_gen, neg_gen2, lclass, None, None]},
                               'LPI': {'[M-H]-': [neg_gen, neg_gen2, lclass, None, None]},
                               'LPS': {'[M-H]-': [neg_gen, neg_gen2, lclass, None, None]},
                               'Cer': {'[M+H]+': [cer_base_spec, cer_fa_spec, 'LCB', None, None]}}
        # Note: add also the length of the list that it should have (georgia: 18.2.2019)
        # Note: in the below section we can add also the info about the link. (georgia: 18.2.2019)
        self.lipid_fa_dct = {'PA': ['PL', ['FA1', 'FA2'], None, None, 2, False],
                             'PC': ['PL', ['FA1', 'FA2'], None, None, 2, True],
                             'PE': ['PL', ['FA1', 'FA2'], None, None, 2, True],
                             'PG': ['PL', ['FA1', 'FA2'], None, None, 2, False],
                             'PI': ['PL', ['FA1', 'FA2'], None, None, 2, False],
                             'PS': ['PL', ['FA1', 'FA2'], None, None, 2, False],
                             'TG': ['TG', ['FA1', 'FA2', 'FA3'], None, None, 3, True],
                             'DG': ['DG', ['FA1', 'FA2'], None, None, 2, True],
                             'LPA': ['LPL', ['FA1'], None, None, 1, False],
                             'LPC': ['LPL', ['FA1'], None, None, 1, True],
                             'LPE': ['LPL', ['FA1'], None, None, 1, True],
                             'LPG': ['LPL', ['FA1'], None, None, 1, False],
                             'LPI': ['LPL', ['FA1'], None, None, 1, False],
                             'LPS': ['LPL', ['FA1'], None, None, 1, False],
                             'CL': ['CL', ['FA1', 'FA2', 'FA3', 'FA4'], None, None, 4, False],
                             'Cer': ['CER', ['FA1'], 'Base', ['M', 'D', 'T'], 2, False]}
        self.lipid_link_dct = {'PA': False}

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

    # Change: move to a futher section by combine it with def get_fa_info and rename it with that name (georgia: 15.1.2019)

    # def calc_fa_all_mz(self, fa_info_dct):
    #
    #     # Change: (georgia_10.1.2019): formation of unecessary function since is used only once
    #     exactmass = self.calc_fa_mass(fa_info_dct)
    #     exit()
    #     fa_info_dct['EXACTMASS'] = exactmass
    #     fa_info_dct['[FA-H]-_MZ'] = round(exactmass - self.elem_dct['H'][0], 6)
    #     fa_info_dct['[FA-H2O-H]-_MZ'] = round(exactmass - 3 * self.elem_dct['H'][0] - self.elem_dct['O'][0], 6)
    #     fa_info_dct['[FA-H2O]_MZ'] = round(exactmass - 2 * self.elem_dct['H'][0] - self.elem_dct['O'][0], 6)
    #     fa_info_dct['[FA-H2O+H]+_MZ'] = round(exactmass - self.elem_dct['H'][0] - self.elem_dct['O'][0], 6)
    #     fa_info_dct['[FA-H+Na]_MZ'] = round(exactmass - self.elem_dct['H'][0] + self.elem_dct['Na'][0], 6)
    #
    #     fa_info_dct['[FA-H]-_ABBR'] = '[{fa}-H]-'.format(fa=fa_info_dct['ABBR'])
    #     fa_info_dct['[FA-H2O-H]-_ABBR'] = '[{fa}-H2O-H]-'.format(fa=fa_info_dct['ABBR'])
    #     fa_info_dct['[FA-H2O]_ABBR'] = '[{fa}-H2O]'.format(fa=fa_info_dct['ABBR'])
    #     fa_info_dct['[FA-H2O+H]+_ABBR'] = '[{fa}-H2O+H]+'.format(fa=fa_info_dct['ABBR'])
    #     fa_info_dct['[FA-H+Na]_ABBR'] = '[{fa}-H+Na]'.format(fa=fa_info_dct['ABBR'])

        # return fa_info_dct
    # Delet: (georgia_10.1.2019)
    # def decode_fa(self, fa_str):
    #
    #     fa_info_lst = []
    #     # Change: (georgia_10.1.2019):
    #     # for _rgx in self.fa_rgx_lst:
    #     # if re.match(_rgx, fa_str):
    #     #     _fa_match = re.match(_rgx, fa_str)
    #     #     fa_info_lst = _fa_match.groups()
    #     #     break
    #
    #     if re.match(self.fa_rgx, fa_str):
    #         _fa_match = re.match(self.fa_rgx, fa_str)
    #         fa_info_lst = _fa_match.groups()
    #
    #     return fa_info_lst

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
        elif fa_info_lst[0] in ['d', 'm', 't']:
            hydroxy_abbr_dct = {'m': 1, 'd': 2, 't': 3}
            fa_info_dct['LINK'] = fa_info_lst[0]
            fa_info_dct['ABBR'] = fa_str
            fa_info_dct['C'] = int(fa_info_lst[1])
            fa_info_dct['H'] = 2 * int(fa_info_lst[1]) - 2 * int(fa_info_lst[3]) + 3
            fa_info_dct['O'] = hydroxy_abbr_dct[fa_info_lst[0]]
            fa_info_dct['N'] = 1
            fa_info_dct['DB'] = fa_info_lst[3]
            fa_info_dct['FORMULA'] = 'C{num_c}H{num_h}O{num_o}N'.format(num_c=fa_info_dct['C'], num_h=fa_info_dct['H'],
                                                                       num_o=fa_info_dct['O'])

        return fa_info_dct

    def get_fa_info(self, fa_str):

        print (fa_str)
        fa_info_dct = self.get_fa_formula(fa_str)
        exactmass = self.calc_fa_mass(fa_info_dct)

        fa_info_dct['EXACTMASS'] = exactmass
        fa_info_dct['[FA-H]-_MZ'] = round(exactmass - self.elem_dct['H'][0], 6)
        fa_info_dct['[FA-H2O-H]-_MZ'] = round(exactmass - 3 * self.elem_dct['H'][0] - self.elem_dct['O'][0], 6)
        fa_info_dct['[FA-H2O]_MZ'] = round(exactmass - 2 * self.elem_dct['H'][0] - self.elem_dct['O'][0], 6)
        fa_info_dct['[FA-H2O+H]+_MZ'] = round(exactmass - self.elem_dct['H'][0] - self.elem_dct['O'][0], 6)
        fa_info_dct['[FA-H+Na]_MZ'] = round(exactmass - self.elem_dct['H'][0] + self.elem_dct['Na'][0], 6)

        fa_info_dct['[FA-H]-_ABBR'] = '[{fa}-H]-'.format(fa=fa_info_dct['ABBR'])
        fa_info_dct['[FA-H2O-H]-_ABBR'] = '[{fa}-H2O-H]-'.format(fa=fa_info_dct['ABBR'])
        fa_info_dct['[FA-H2O]_ABBR'] = '[{fa}-H2O]'.format(fa=fa_info_dct['ABBR'])
        fa_info_dct['[FA-H2O+H]+_ABBR'] = '[{fa}-H2O+H]+'.format(fa=fa_info_dct['ABBR'])
        fa_info_dct['[FA-H+Na]_ABBR'] = '[{fa}-H+Na]'.format(fa=fa_info_dct['ABBR'])

        return fa_info_dct


    # Compine with the calc_fa_mass as above
    # def get_fa_info(self, fa_str):
    #     print (fa_str)
    #     fa_info_dct = self.get_fa_formula(fa_str)
    #     fa_info_dct = self.calc_fa_all_mz(fa_info_dct)
    #
    #     return fa_info_dct

    # Note: More general way to define the fragments. Also define specific and not. (georgia: 23.1.2019)
    #Change: Values for HIGH and LOW can be removed. Check if are need or not (georgia: 15.2.2019)
    # TODO: check if they cause an issue or not. (low_v & high_v) (georgia: 15.2.2019)
    def frag_definition(self, lclass, mode, fa_info , ms2_ppm = 20):
        print ('lets see what we are going to do here')
        exactmass = self.calc_fa_mass(fa_info)
        fa_info['EXACTMASS'] = exactmass
        low_v = ''
        high_v = ''
        if fa_info['LINK'] in ['d', 'm', 't']:
            print ('welll')
            fragment_list = self.lipid_mode_dct[lclass][mode][0]

            for frag in fragment_list:
                fa_info['{f}_MZ'.format(f=frag)]= round(exactmass + self.loss_dct[frag], 6)
                fa_info['{f}_ABBR'.format(f=frag)] = frag.replace('LCB', 'LCB({l})'.format(l=fa_info['ABBR']))
                low_v = round(ppm_window_para(fa_info['{f}_MZ'.format(f=frag)], ms2_ppm * -1), 6)
                high_v = round(ppm_window_para(fa_info['{f}_MZ'.format(f=frag)], ms2_ppm), 6)
                fa_info['{f}_Q'.format(f=frag)] = (low_v.astype(str) + ' <= mz <= ' + high_v.astype(str))


        elif fa_info['LINK'] in ['FA', 'O-', 'P-']:
            if lclass == 'Cer':
                fragment_list = self.lipid_mode_dct[lclass][mode][1]
                for frag in fragment_list:
                    fa_info['{f}_MZ'.format(f=frag)] = round(exactmass + self.loss_dct[frag], 6)
                    fa_info['{f}_ABBR'.format(f=frag)] = frag.replace('FA', fa_info['ABBR'])
                    low_v = round(ppm_window_para(fa_info['{f}_MZ'.format(f=frag)],
                                                                           ms2_ppm * -1), 6)
                    high_v = round(ppm_window_para(fa_info['{f}_MZ'.format(f=frag)], ms2_ppm), 6)
                    fa_info['{f}_Q'.format(f=frag)] = (low_v.astype(str)) + ' <= mz <= ' + (high_v.astype(str))
            else:
                fragment_list = self.lipid_mode_dct[lclass][mode][0]
                for frag in fragment_list:
                    fa_info['{f}_MZ'.format(f=frag)] = round(exactmass + self.loss_dct[frag], 6)
                    fa_info['{f}_ABBR'.format(f=frag)] = frag.replace('FA', fa_info['ABBR'])
                    low_v = round(ppm_window_para(fa_info['{f}_MZ'.format(f=frag)], ms2_ppm * -1), 6)
                    high_v = round(ppm_window_para(fa_info['{f}_MZ'.format(f=frag)], ms2_ppm), 6)
                    fa_info['{f}_Q'.format(f=frag)] = (low_v.astype(str) + ' <= mz <= ' + high_v.astype(str))

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
    def get_fa_info_geo(self, fa_str, lipid_class, ion_mode, ms2_ppm = 20):

        print (fa_str)
        fa_info_dct = self.get_fa_formula(fa_str)
        fa_info_dct = self.frag_definition(lipid_class, ion_mode, fa_info_dct)


        return fa_info_dct


if __name__ == '__main__':

    abbr_decoder = NameParserFA()
    lipidC = 'Cer'
    ionM = '[M+H]+'
    abbr_lst = ['d18:0', 'FA8:0',  'FA18:2', 'O-16:0', 'P-18:0']
    for abbr in abbr_lst:
        x = abbr_decoder.get_fa_info_geo(abbr, lipidC, ionM)
        print(x)
