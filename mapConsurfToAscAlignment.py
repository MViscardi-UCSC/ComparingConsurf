"""
mapConSurfOnToAlignment.py
Marcus Viscardi  May 16, 2020

Taking parsing scripts from parseASC.py and parseConSurfGrades.py

Using these to map the consurf data onto the alignments
"""
from time import sleep
from typing import Dict, List

import pandas as pd

pd.set_option('display.max_rows', 500)
pd.set_option('display.max_columns', 500)
pd.set_option('display.width', 1000)

DF_DICT = Dict[str, pd.DataFrame]
STR_DICT = Dict[str, str]


def parse_grades(file: str, header_lines: int, footer: int):
    dataframe = pd.read_fwf(open(file, 'r'), skiprows=header_lines, skipfooter=footer)
    pattern = r"(\d+:[A-Z])"
    dataframe['ID'] = dataframe['3LATOM'].str.extract(pat=pattern)
    dataframe['ID'] = ':' + dataframe['ID'].str.replace(":", ".")
    output = dataframe[['ID', 'SEQ', 'SCORE']]
    print(output.head(10))
    return output


def double_grades_wrapper(files: STR_DICT, header: int, footer: int) -> DF_DICT:
    df_dict = {}
    for name, file in files.items():
        df = parse_grades(file, header, footer)
        df_dict[name] = df
        print(f"{name}: {df.shape[0]} graded residues")
    return df_dict


def parse_asc(file: str, split_line_num: int) -> DF_DICT:
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
        df_dict[name] = df[['POS', 'ID']]
        print(f"{name}: {df.shape[0]} aligned residues")
    return df_dict


def map_together(alignment_dict: DF_DICT, grades_dict: DF_DICT):
    print('\nMerging  dataframes:\n')
    final_df_dict = {}
    for name, df in grades_dict.items():
        try:
            # This merge brings together the grades data, containing the consurf scores with the alignment data.
            #   They are matched together based on the unique residue IDs
            final_df = df.merge(alignment_dict[name], on="ID", how='left')
            final_df_dict[name] = final_df
        except KeyError as ke:
            print(f"ERROR: Mismatch in dataframe dictionary keys. {ke} "
                  f"not found in keys of alignment dictionary: {list(alignment_dict.keys())}")
            sleep(0.2)
            raise ke
    key_list = [k for k in grades_dict.keys()]
    HWX, HWY = key_list
    # This second merge matches the two dataframes based on the POS positional values, which have been carried over
    #   from the original alignment files. This matching allows for us to find the difference in Consurf scores!
    merge_df = final_df_dict[HWX].merge(final_df_dict[HWY], on='POS', how='left',
                                        suffixes=(f'_{HWX}', f'_{HWY}'))
    merge_df['SCORE_DIFF'] = pd.Series(merge_df[f'SCORE_{HWX}'] - merge_df[f'SCORE_{HWY}']).abs()
    output_df_dict = {}
    for name in key_list:
        output_df_dict[name] = merge_df
        print(output_df_dict[name].head(10))
    return output_df_dict


def output_attribute_files(df_dict: DF_DICT, suffix: str = '', attribute_name: str = 'consurfScoreDiff') -> None:
    key_list = [k for k in df_dict.keys()]
    header_lines = [f"# Comparing consurf scores of {key_list[0]} and {key_list[1]}",
                    f"attribute: {attribute_name}",
                    "match mode: 1-to-1",
                    "recipient: residues"]
    for name, df in df_dict.items():
        # Output regular Consurf results:
        output_file = f'{name}_ConSurfScore.txt'
        attribute_lines = []
        for line in df.to_string(columns=[f'ID_{name}', f'SCORE_{name}'], index=False).split('\n'):
            if not line.strip().startswith(('ID', 'NaN')):
                id_att_list = line.strip().split()
                line = f"\t{id_att_list[0]}\t{id_att_list[1]}"
                attribute_lines.append(line)
        # print(name, attribute_lines, sep='\n')
        all_lines = [f"# Consurf scores of {name}",
                     f"attribute: consurfScore",
                     "match mode: 1-to-1",
                     "recipient: residues"] + attribute_lines
        all_lines = [f'{l}\n' for l in all_lines]
        with open(output_file, 'w') as file:
            file.writelines(all_lines)

        # Regular diff output
        attribute_lines = []
        output_file = f'{name}_ConSurfScoreDiff{suffix}.txt'
        for line in df.to_string(columns=[f'ID_{name}', f'SCORE_DIFF'], index=False).split('\n'):
            if 'NaN' not in line.strip() and 'ID' not in line.strip():
                id_att_list = line.strip().split()
                line = f"\t{id_att_list[0]}\t{id_att_list[1]}"
                attribute_lines.append(line)
        # print(name, attribute_lines, sep='\n')
        all_lines = header_lines + attribute_lines
        all_lines = [f'{l}\n' for l in all_lines]
        with open(output_file, 'w') as file:
            file.writelines(all_lines)


def main(grades_dict: STR_DICT, asc_path_and_breakline: List[str and int], suffix: str = ''):
    alignment_dict = parse_asc(asc_path_and_breakline[0], asc_path_and_breakline[1])
    grades_dict = double_grades_wrapper(grades_dict, 13, 4)
    output_df = map_together(alignment_dict, grades_dict)
    output_attribute_files(output_df, suffix=suffix)


if __name__ == '__main__':
    used_grades_dict = {"2HWX": r"T:\Chrome Downloads\SMG5and6\SMG6_(2HWX)\All_Consurf_Outputs\consurf.grades",
                        "2HWY": r"T:\Chrome Downloads\SMG5and6\SMG5_(2HWY)\All_Consurf_Outputs\consurf.grades", }
    asc_path = [r"T:\Chrome Downloads\SMG5and6\2HWY&2HWX\Realigned Sequences.asc", 167]
    asc_path_no_realign = [r"T:\Chrome Downloads\SMG5and6\2HWY&2HWX\2HWXand2HWY_notRealigned.asc", 118]
    main(used_grades_dict, asc_path, suffix='_realigned')
    main(used_grades_dict, asc_path_no_realign, suffix='_noRealign')
