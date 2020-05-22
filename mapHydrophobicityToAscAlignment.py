"""
mapHydrophobicityToASCAlignment.py
Marcus Viscardi  May 21, 2020

"""
import sys
from typing import Dict
import pandas as pd

pd.set_option('display.max_rows', 500)
pd.set_option('display.max_columns', 500)
pd.set_option('display.width', 1000)

DF_DICT = Dict[str, pd.DataFrame]


def parse_three_letter_asc(file: str, split_line_num: int) -> DF_DICT:
    with open(file, "r") as f:
        df1 = pd.read_table(f, nrows=split_line_num - 2)
    with open(file, "r") as f:
        df2 = pd.read_table(f, skiprows=split_line_num - 1)
    both = [df1, df2]
    df_dict = {}
    for df in both:
        column = df.columns[0]
        name = column[:4]
        # It is important that I am keeping this POS positional data, not the position from the .grades files!
        df[['POS', 'ID']] = pd.DataFrame(df.apply(lambda x: x[column].split(' - '), axis=1).tolist(), index=df.index)
        df = df.astype({column: 'object',
                        'POS': 'int64',
                        'ID': 'object'})
        df[['RES', 'ID']] = pd.DataFrame(df.apply(lambda x: x['ID'].split(), axis=1).tolist(), index=df.index)
        df['ID'] = ':' + df['ID']
        df_dict[name] = df[['POS', 'RES', 'ID']]
        print(f"{name}: {df.shape[0]} aligned residues")
    return df_dict


def add_hydrophobicity(df_dict: DF_DICT, scale: str = "kd"):
    """
    Hydrophobicity scales: (https://www.cgl.ucsf.edu/chimera/docs/UsersGuide/midas/hydrophob.html)
        RES KD  WW      HH      MF      TT
        Ile	4.5	0.31	-0.60	-1.56	1.97
        Val	4.2	-0.07	-0.31	-0.78	1.46
        Leu	3.8	0.56	-0.55	-1.81	1.82
        Phe	2.8	1.13	-0.32	-2.20	1.98
        Cys	2.5	0.24	-0.13	0.49	-0.30
        Met	1.9	0.23	-0.10	-0.76	1.40
        Ala	1.8	-0.17	0.11	0.0	0.38
        Gly	-0.4	-0.01	0.74	1.72	-0.19
        Thr	-0.7	-0.14	0.52	1.78	-0.32
        Ser	-0.8	-0.13	0.84	1.83	-0.53
        Trp	-0.9	1.85	0.30	-0.38	1.53
        Tyr	-1.3	0.94	0.68	-1.09	0.49
        Pro	-1.6	-0.45	2.23	-1.52	-1.44
        His	-3.2	-0.96	2.06	4.76	-1.44
        Glu	-3.5	-2.02	2.68	1.64	-2.90
        Gln	-3.5	-0.58	2.36	3.01	-1.84
        Asp	-3.5	-1.23	3.49	2.95	-3.27
        Asn	-3.5	-0.42	2.05	3.47	-1.62
        Lys	-3.9	-0.99	2.71	5.39	-3.46
        Arg	-4.5	-0.81	2.58	3.71	-2.57
    """
    scales = {'kd': {'ile': 4.5,
                     'val': 4.2,
                     'leu': 3.8,
                     'phe': 2.8,
                     'cys': 2.5,
                     'met': 1.9,
                     'ala': 1.8,
                     'gly': -0.4,
                     'thr': -0.7,
                     'ser': -0.8,
                     'trp': -0.9,
                     'tyr': -1.3,
                     'pro': -1.6,
                     'his': -3.2,
                     'glu': -3.5,
                     'gln': -3.5,
                     'asp': -3.5,
                     'asn': -3.5,
                     'lys': -3.9,
                     'arg': -4.5,
                     },
              'ww': {'ile': 0.31,
                     'val': -0.07,
                     'leu': 0.56,
                     'phe': 1.13,
                     'cys': 0.24,
                     'met': 0.23,
                     'ala': -0.17,
                     'gly': -0.01,
                     'thr': -0.14,
                     'ser': -0.13,
                     'trp': 1.85,
                     'tyr': 0.94,
                     'pro': -0.45,
                     'his': -0.96,
                     'glu': -2.02,
                     'gln': -0.58,
                     'asp': -1.23,
                     'asn': -0.42,
                     'lys': -0.99,
                     'arg': -0.81,
                     },
              'hh': {'ile': -0.60,
                     'val': -0.31,
                     'leu': -0.55,
                     'phe': -0.32,
                     'cys': -0.13,
                     'met': -0.10,
                     'ala': 0.11,
                     'gly': 0.74,
                     'thr': 0.52,
                     'ser': 0.84,
                     'trp': 0.30,
                     'tyr': 0.68,
                     'pro': 2.23,
                     'his': 2.06,
                     'glu': 2.68,
                     'gln': 2.36,
                     'asp': 3.49,
                     'asn': 2.05,
                     'lys': 2.71,
                     'arg': 2.58,
                     },
              }
    if scale.lower() == "kd":  # TODO: This would be better to have each be a key to a nested dict
        scale = {'ile': 4.5,
                 'val': 4.2,
                 'leu': 3.8,
                 'phe': 2.8,
                 'cys': 2.5,
                 'met': 1.9,
                 'ala': 1.8,
                 'gly': -0.4,
                 'thr': -0.7,
                 'ser': -0.8,
                 'trp': -0.9,
                 'tyr': -1.3,
                 'pro': -1.6,
                 'his': -3.2,
                 'glu': -3.5,
                 'gln': -3.5,
                 'asp': -3.5,
                 'asn': -3.5,
                 'lys': -3.9,
                 'arg': -4.5,
                 }
    else:  # Then could catch this as an IndexError
        print(f"Scale: '{scale}' not currently accepted, sorry")
        sys.exit()
    # First add the hydrophobicity scores on from the scale dictionary:
    for name, df in df_dict.items():
        df['HYDRO'] = pd.Series(df['RES'].str.lower()).map(scale)
    # Then merge the two dicts together based on the alignment POS
    key_list = [k for k in df_dict.keys()]
    hwy, hwx = key_list
    print(f"\nKeys:\n\tHWX: '{hwx}'\n\tHWY: '{hwy}'")
    merge_df = df_dict[hwx].merge(df_dict[hwy], on='POS', how='outer', suffixes=(f'_{hwx}', f'_{hwy}'))
    merge_df['HYDRO_DIFF'] = pd.Series(merge_df[f'HYDRO_{hwx}'] - merge_df[f'HYDRO_{hwy}']).abs()
    print(merge_df[['POS', f'ID_{hwx}', f'ID_{hwy}', 'HYDRO_DIFF']])
    output_df_dict = {}
    for name in key_list:
        output_df_dict[name] = merge_df[[f'ID_{name}', f'HYDRO_{name}', 'HYDRO_DIFF']]
        print(output_df_dict[name].head(10))
    return output_df_dict


def output_attribute_files(df_dict: DF_DICT, suffix: str = '',
                           scale: str = 'kd') -> None:
    key_list = [k for k in df_dict.keys()]
    header_lines = [f"# Comparing {scale}Hydrophobicity values of {key_list[0]} and {key_list[1]}",
                    f"attribute: {scale}_Hydro_Diff",
                    "match mode: 1-to-1",
                    "recipient: residues"]
    for name, df in df_dict.items():
        # Regular diff output
        attribute_lines = []
        output_file = f'{name}_{scale}HydrophobicityDiff{suffix}.txt'
        for line in df.to_string(columns=[f'ID_{name}', f'HYDRO_DIFF'], index=False).split('\n'):
            if 'NaN' not in line.strip() and 'ID' not in line.strip():
                id_att_list = line.strip().split()
                line = f"\t{id_att_list[0]}\t{id_att_list[1]}"
                attribute_lines.append(line)
        # print(name, attribute_lines, sep='\n')
        all_lines = header_lines + attribute_lines
        all_lines = [f'{l}\n' for l in all_lines]
        open(output_file, 'w').writelines(all_lines)


if __name__ == '__main__':
    file = r"Realigned_Sequences_w3letter.asc"
    dataframe_dict = parse_three_letter_asc(file, 118)
    scale = 'kd'
    diff_df_dict = add_hydrophobicity(dataframe_dict, scale=scale)
    output_attribute_files(diff_df_dict, scale=scale)
