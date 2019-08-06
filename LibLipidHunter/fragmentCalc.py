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

try:
    from LibLipidHunter.LipidNomenclature import NameParserFA
    from LibLipidHunter.AbbrElemCalc import ElemCalc
    from LibLipidHunter.ParallelFunc import ppm_window_para
except ImportError:  # for python 2.7.14
    from LipidNomenclature import NameParserFA
    from AbbrElemCalc import ElemCalc
    from ParallelFunc import ppm_window_para

class allFrag:
    def __init__(self, lclass):
        # ElemCalc().lipid_hg_elem_dct
        # ElemCalc().glycerol_bone_elem_dct
        # Note: second line is to make it more genreal
        elemen_mass = NameParserFA().elem_dct
        neg_gen = ['[FA-H]-', '[FA-H2O]', '[FA-H2O-H]-']
        neg_gen2 = ['[x-H]-', '[x-H2O-H]-']
        pos_gen = ['[FA-H2O+H]+', '[FA-H2O]']
        pos_gen2 = ['[x-H2O+H]+', '[x-H2O]']
        sodium_spec = ['[FA-H+Na]']
        sodium_spec2 = ['[x-H+Na]']
        cer_base_spec = ['[LCB-H2O+H]+', '[LCB-2xH2O+H]+', '[LCB-H2O-CH2O+H]+']
        cer_base_spec2 = ['[x-H2O+H]+', '[x-2xH2O+H]+', '[x-H2O-CH2O+H]+']
        cer_fa_spec = ['[FA-OH+NH3]+', '[FA-OH+C2H3N]+']
        cer_fa_spec2 = ['[x-OH+NH3]+', '[x-OH+C2H3N]+']

        if lclass == 'PC':
            formicacid_loss = -(elemen_mass['C'][0] + elemen_mass['H'][0]*3)
        else:
            formicacid_loss = 0

        self.loss_dct = {'[FA-H]-': -elemen_mass['H'][0],
                         '[FA-H2O]': -(elemen_mass['H'][0]*2 + elemen_mass['O'][0]),
                         '[FA-H2O-H]-': -(elemen_mass['H'][0]*3 + elemen_mass['O'][0]),
                         '[FA-H2O+H]+': -(elemen_mass['H'][0] + elemen_mass['O'][0]),
                         '[FA-H+Na]': elemen_mass['Na'][0] - elemen_mass['H'][0],
                         '[LCB-H2O+H]+': -(elemen_mass['H'][0] + elemen_mass['O'][0]),
                         '[LCB-2xH2O+H]+': -(elemen_mass['H'][0]*3 + elemen_mass['O'][0]*2),
                         '[LCB-H2O-CH2O+H]+': -(elemen_mass['H'][0]*3 + elemen_mass['O'][0]*2 + elemen_mass['C'][0]),
                         '[FA-OH+NH3]+': elemen_mass['N'][0] + elemen_mass['H'][0]*2 - elemen_mass['O'][0],
                         '[FA-OH+C2H3N]+': (elemen_mass['C'][0] + elemen_mass['N'][0] + elemen_mass['H'][0]*2 -
                                            elemen_mass['O'][0]),
                         '[x-H]-': NameParserFA().calc_fa_mass(ElemCalc().lipid_hg_elem_dct[lclass])
                                   + NameParserFA().calc_fa_mass(ElemCalc().glycerol_bone_elem_dct)
                                   + elemen_mass['H'][0] + elemen_mass['O'][0] + formicacid_loss,
                         '[x-H2O-H]-': NameParserFA().calc_fa_mass(ElemCalc().lipid_hg_elem_dct[lclass])
                                   + NameParserFA().calc_fa_mass(ElemCalc().glycerol_bone_elem_dct)
                                   - elemen_mass['H'][0] + formicacid_loss}

        self.lipid_mode_dct = {'PA': {'[M-H]-': [neg_gen, neg_gen2, 'LPA']},
                               'PC': {'[M+HCOO]-': [neg_gen, neg_gen2, 'LPC'], '[M+CH3COO]-': [neg_gen, neg_gen2, 'LPC']},
                               'PE': {'[M-H]-': [neg_gen, neg_gen2, 'LPE']},
                               'PG': {'[M-H]-': [neg_gen, neg_gen2, 'LPG']},
                               'PI': {'[M-H]-': [neg_gen, neg_gen2, 'LPI']},
                               'PS': {'[M-H]-': [neg_gen, neg_gen2, 'LPS']},
                               'TG': {'[M+NH4]+': [pos_gen, [], ''], '[M+H]+': [pos_gen, [], ''], '[M+Na]+': [(pos_gen + sodium_spec), [], '']},
                               'DG': {'[M+NH4]+': [pos_gen, [], '']},
                               'LPA': {'[M-H]-': [neg_gen, neg_gen2, lclass]},
                               'LPC': {'[M+HCOO]-': [neg_gen, neg_gen2, lclass], '[M+CH3COO]-': [neg_gen, neg_gen2, lclass]},
                               'LPE': {'[M-H]-': [neg_gen, neg_gen2, lclass]},
                               'LPG': {'[M-H]-': [neg_gen, neg_gen2, lclass]},
                               'LPI': {'[M-H]-': [neg_gen, neg_gen2, lclass]},
                               'LPS': {'[M-H]-': [neg_gen, neg_gen2, lclass]},
                               'Cer': {'[M+H]+': [cer_base_spec, cer_fa_spec, 'LCB']}}
        #print (self.lipid_mode_dct)

    def frag_definition(self, lclass, mode, fa_info , ms2_ppm = 20):
        print ('lets see what we are going to do here')
        exactmass = NameParserFA().calc_fa_mass(fa_info)
        fa_info['EXACTMASS'] = exactmass
        if fa_info['LINK'] in ['d', 'm', 't']:
            print ('welll')
            fragment_list = self.lipid_mode_dct[lclass][mode][0]
            for frag in fragment_list:
                fa_info['{f}_MZ'.format(f=frag)]= exactmass + self.loss_dct[frag]
                fa_info['{f}_ABBR'.format(f=frag)] = frag.replace('LCB', 'LCB({l})'.format(l=fa_info['ABBR']))

        elif fa_info['LINK'] in ['FA', 'O-', 'P-']:
            if lclass == 'Cer':
                fragment_list = self.lipid_mode_dct[lclass][mode][1]
                for frag in fragment_list:
                    fa_info['{f}_MZ'.format(f=frag)] = exactmass + self.loss_dct[frag]
                    fa_info['{f}_ABBR'.format(f=frag)] = frag.replace('FA', fa_info['ABBR'])
            else:
                fragment_list = self.lipid_mode_dct[lclass][mode][0]
                for frag in fragment_list:
                    fa_info['{f}_MZ'.format(f=frag)] = exactmass + self.loss_dct[frag]
                    fa_info['{f}_ABBR'.format(f=frag)] = frag.replace('FA', fa_info['ABBR'])
                    fa_info['{f}_MZ_LOW'.format(f=frag)] = ppm_window_para(fa_info['{f}_MZ'.format(f=frag)], ms2_ppm * -1)
                    fa_info['{f}_MZ_HIGH'.format(f=frag)] = ppm_window_para(fa_info['{f}_MZ'.format(f=frag)], ms2_ppm)
                    fa_info['{f}_Q'.format(f=frag)] = (fa_info['{f}_MZ_LOW'.format(f=frag)].astype(str)) + ' <= mz <= ' + (fa_info['{f}_MZ_HIGH'.format(f=frag)].astype(str))

                fragment_list2 = self.lipid_mode_dct[lclass][mode][1]
                fragm2_name = "{l}({f})".format(l=self.lipid_mode_dct[lclass][mode][2], f=fa_info['ABBR'].replace('FA', ''))

                for frag in fragment_list2:
                    fragm_name = frag.replace('x', self.lipid_mode_dct[lclass][mode][2])
                    fa_info['{f}_MZ'.format(f=fragm_name)] = exactmass + self.loss_dct[frag]
                    fa_info['{f}_ABBR'.format(f=fragm_name)] = frag.replace('x', fragm2_name)
                    fa_info['{f}_MZ_LOW'.format(f=fragm_name)] = ppm_window_para(fa_info['{f}_MZ'.format(f=fragm_name)],
                                                                           ms2_ppm * -1)
                    fa_info['{f}_MZ_HIGH'.format(f=fragm_name)] = ppm_window_para(fa_info['{f}_MZ'.format(f=fragm_name)], ms2_ppm)
                    fa_info['{f}_Q'.format(f=fragm_name)] = (fa_info['{f}_MZ_LOW'.format(f=fragm_name)].astype(
                        str)) + ' <= mz <= ' + (fa_info['{f}_MZ_HIGH'.format(f=fragm_name)].astype(str))


        return (fa_info)

if __name__ == '__main__':


    info_class = 'PC'
    info_mode = '[M+HCOO]-'
    abbr_decoder = allFrag(info_class)
    # info_dct={'LINK': 'd', 'ABBR': 'd6:0', 'C': 6, 'H': 15, 'O': 2, 'N': 1, 'DB': '0', 'FORMULA': 'C6H15O2N'}
    # info_dct = {'LINK': 'd', 'ABBR': 'd18:0', 'C': 18, 'H': 39, 'O': 2, 'N': 1, 'DB': '0', 'FORMULA': 'C18H39O2N'}
    info_dct = {'ABBR': 'FA18:0', 'LINK': 'FA', 'C': 18, 'H': 36, 'O': 2, 'DB': 0, 'FORMULA': 'C18H36O2'}
    info_dct = abbr_decoder.frag_definition(info_class, info_mode, info_dct)
    print (info_dct)