"""
ConsurfScatterplot.py
Marcus Viscardi     May 23, 2020

Using example from:
https://medium.com/@gorjanz/data-analysis-in-python-interactive-scatterplot-with-matplotlib-6bb8ad2f1f18

Going to use some of the annotation parsing from mapConsurf... and mapHydro...

Goal is to make a plot with the conservation in 2HWY on one axis, and 2HWX on the other. Each point will be an aligned
    pair of residues. By clicking on the residues in the plot it will annotate the two relevant residues
"""
import sys

import pandas as pd
from typing import Dict, List
import matplotlib.pyplot as plt
import numpy as np

from mapHydrophobicityToAscAlignment import add_hydrophobicity

pd.set_option('display.max_rows', 500)
pd.set_option('display.max_columns', 500)
pd.set_option('display.width', 1000)

plt.style.use('ggplot')

DF_DICT = Dict[str, pd.DataFrame]
STR_DICT = Dict[str, str]


def example():
    global generated_labels, instances_colors, axis_values_x, axis_values_y, ax, draw_scatterplot, annotate
    # import the random module since we will use it to generate the data
    import random as rnd
    # import the main drawing library
    import matplotlib.pyplot as plt
    from matplotlib.widgets import Button
    from matplotlib.text import Annotation
    # import the seaborn module which is based on matplotlib to make our visualization more presentable
    import seaborn as sns
    # set the default style
    sns.set()
    # define two colors, just to enrich the example
    labels_color_map = {0: '#20b2aa', 1: '#ff7373'}
    # set the examples count
    no_examples = 50
    # generate the data needed for the scatterplot
    generated_data = [(x, rnd.randint(0, no_examples)) for x in range(0, no_examples)]
    generated_labels = ["Label for instance #{0}".format(i) for i in range(0, no_examples)]
    print("now visualizing scatterlplot...")
    # add the values one by one to the scatterplot
    instances_colors = []
    axis_values_x = []
    axis_values_y = []
    for index, instance in enumerate(generated_data):
        # print instance, index, labels[index]
        coordinate_x, coordinate_y = instance
        color = labels_color_map[index % 2]

        instances_colors.append(color)
        axis_values_x.append(coordinate_x)
        axis_values_y.append(coordinate_y)
    # draw a scatter-plot of the generated values
    fig = plt.figure(figsize=(20, 16))
    ax = plt.subplot()

    # extract the scatterplot drawing in a separate function so we ca re-use the code
    def draw_scatterplot():
        ax.scatter(
            axis_values_x,
            axis_values_y,
            c=instances_colors,
            picker=True
        )

    # draw the initial scatterplot
    draw_scatterplot()

    # create and add an annotation object (a text label)
    def annotate(axis, text, x, y):
        text_annotation = Annotation(text, xy=(x, y), xycoords='data')
        axis.add_artist(text_annotation)

    # define the behaviour -> what happens when you pick a dot on the scatterplot by clicking close to it
    def onpick(event):
        # step 1: take the index of the dot which was picked
        ind = event.ind

        # step 2: save the actual coordinates of the click, so we can position the text label properly
        label_pos_x = event.mouseevent.xdata
        label_pos_y = event.mouseevent.ydata

        # just in case two dots are very close, this offset will help the labels not appear one on top of each other
        offset = 0

        # if the dots are to close one to another, a list of dots clicked is returned by the matplotlib library
        for i in ind:
            # step 3: take the label for the corresponding instance of the data
            label = generated_labels[i]

            # step 4: log it for debugging purposes
            print("index", i, label)

            # step 5: create and add the text annotation to the scatterplot
            annotate(
                ax,
                label,
                label_pos_x + offset,
                label_pos_y + offset
            )

            # step 6: force re-draw
            ax.figure.canvas.draw_idle()

            # alter the offset just in case there are more than one dots affected by the click
            offset += 0.01

    # connect the click handler function to the scatterplot
    fig.canvas.mpl_connect('pick_event', onpick)
    # create the "clear all" button, and place it somewhere on the screen
    ax_clear_all = plt.axes([0.0, 0.0, 0.1, 0.05])
    button_clear_all = Button(ax_clear_all, 'Clear all')

    # define the "clear all" behaviour
    def onclick(event):
        # step 1: we clear all artist object of the scatter plot
        ax.cla()

        # step 2: we re-populate the scatterplot only with the dots not the labels
        draw_scatterplot()

        # step 3: we force re-draw
        ax.figure.canvas.draw_idle()

    # link the event handler function to the click event on the button
    button_clear_all.on_clicked(onclick)
    # initial drawing of the scatterplot
    plt.plot()
    print("scatterplot done")
    # present the scatterplot
    plt.show()


def my_denovo_hydro(asc_file: str, split_line: int, scale: str) -> pd.DataFrame:
    def parse_three_letter_asc(file: str, split_line_num: int) -> DF_DICT:
        with open(file, "r") as f:
            df1 = pd.read_csv(f, sep='\t', nrows=split_line_num - 2)
        with open(file, "r") as f:
            df2 = pd.read_csv(f, sep='\t', skiprows=split_line_num - 1)
        both = [df1, df2]
        df_dict = {}
        for df in both:
            column = df.columns[0]
            name = column[:4]
            # It is important that I am keeping this POS positional data, not the position from the .grades files!
            df[['POS', 'ID']] = pd.DataFrame(df.apply(lambda x: x[column].split(' - '), axis=1).tolist(),
                                             index=df.index)
            df = df.astype({column: 'object',
                            'POS': 'int64',
                            'ID': 'object'})
            df[['RES', 'ID']] = pd.DataFrame(df.apply(lambda x: x['ID'].split(), axis=1).tolist(), index=df.index)
            df['ID'] = ':' + df['ID']
            df_dict[name] = df[['POS', 'RES', 'ID']]
            print(f"{name}: {df.shape[0]} aligned residues")
        return df_dict

    def add_hydrophobicity(df_dict: DF_DICT, scale: str = "kd") -> DF_DICT:
        # Hydrophobicity scales: (https://www.cgl.ucsf.edu/chimera/docs/UsersGuide/midas/hydrophob.html)
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
                  'mf': {'ile': -1.56,
                         'val': -0.78,
                         'leu': -1.81,
                         'phe': -2.20,
                         'cys': 0.49,
                         'met': -0.76,
                         'ala': 0.0,
                         'gly': 1.72,
                         'thr': 1.78,
                         'ser': 1.83,
                         'trp': -0.38,
                         'tyr': -1.09,
                         'pro': -1.52,
                         'his': 4.76,
                         'glu': 1.64,
                         'gln': 3.01,
                         'asp': 2.95,
                         'asn': 3.47,
                         'lys': 5.39,
                         'arg': 3.71,
                         },
                  'tt': {'ile': 1.97,
                         'val': 1.46,
                         'leu': 1.82,
                         'phe': 1.98,
                         'cys': -0.30,
                         'met': 1.40,
                         'ala': 0.38,
                         'gly': -0.19,
                         'thr': -0.32,
                         'ser': -0.53,
                         'trp': 1.53,
                         'tyr': 0.49,
                         'pro': -1.44,
                         'his': -1.44,
                         'glu': -2.90,
                         'gln': -1.84,
                         'asp': -3.27,
                         'asn': -1.62,
                         'lys': -3.46,
                         'arg': -2.57,
                         },
                  }
        if scale.lower() in scales.keys():
            scale = scales[scale]
        else:
            print(f"Scale: '{scale}' not accepted, sorry?")
            sys.exit()
        # Add the hydrophobicity scores on from the scale dictionary:
        for name, df in df_dict.items():
            df['HYDRO'] = pd.Series(df['RES'].str.lower()).map(scale)
        return df_dict

    def merge_df_dict(df_dict: DF_DICT) -> pd.DataFrame:
        key_list = [k for k in df_dict.keys()]
        hwy, hwx = key_list
        print(f"\nKeys:\n\tHWX: '{hwx}'\n\tHWY: '{hwy}'")
        merge_df = df_dict[hwx].merge(df_dict[hwy], on='POS', how='outer', suffixes=(f'_{hwx}', f'_{hwy}'))
        merge_df['HYDRO_DIFF'] = pd.Series(merge_df[f'HYDRO_{hwx}'] - merge_df[f'HYDRO_{hwy}']).abs()
        merge_df = merge_df.dropna()
        merge_df = merge_df.set_index('POS')
        print(merge_df)
        return merge_df

    def main_df_generation(asc_file: str, split_line: int, scale: str) -> pd.DataFrame:
        dataframe_dict = parse_three_letter_asc(asc_file, split_line)
        hydro_df_dict = add_hydrophobicity(dataframe_dict, scale=scale)
        merge_df = merge_df_dict(hydro_df_dict)
        return merge_df

    dataframe = main_df_generation(asc_file=asc_file, split_line=split_line, scale=scale)
    print(dataframe.columns)
    df = dataframe
    fig, ax = plt.subplots(facecolor='w')
    x, y = '2HWX', '2HWY'
    print('\n\n\n\n', df)
    print(df['HYDRO_DIFF'].max())
    ax.scatter(df[f'HYDRO_{x}'], df[f'HYDRO_{y}'],
               c=df['HYDRO_DIFF'],
               cmap='Spectral',
               alpha=0.75,
               )
    # TODO: Add annotations! Maybe only add them to the ones above some threshold
    # line = np.linspace(-5, 5, 11)
    # plt.plot(line, line, 'k-', alpha=0.25)
    ax.set_xlabel(f'Hydrophobicity of residues in {x} on kd scale')
    ax.set_ylabel(f'Hydrophobicity of residues in {y} on kd scale')
    plt.show()
    return dataframe


def my_denovo_consurf(grades_dict: STR_DICT, asc_file: str, asc_split_line: int, suffix: str = ''):
    def parse_three_letter_asc(file: str, split_line_num: int) -> DF_DICT:
        with open(file, "r") as f:
            df1 = pd.read_csv(f, sep='\t', nrows=split_line_num - 2)
        with open(file, "r") as f:
            df2 = pd.read_csv(f, sep='\t', skiprows=split_line_num - 1)
        both = [df1, df2]
        df_dict = {}
        for df in both:
            column = df.columns[0]
            name = column[:4]
            # It is important that I am keeping this POS positional data, not the position from the .grades files!
            df[['POS', 'ID']] = pd.DataFrame(df.apply(lambda x: x[column].split(' - '), axis=1).tolist(),
                                             index=df.index)
            df = df.astype({column: 'object',
                            'POS': 'int64',
                            'ID': 'object'})
            df[['RES', 'ID']] = pd.DataFrame(df.apply(lambda x: x['ID'].split(), axis=1).tolist(), index=df.index)
            df['ID'] = ':' + df['ID']
            df_dict[name] = df[['POS', 'RES', 'ID']]
            print(f"{name}: {df.shape[0]} aligned residues")
        return df_dict

    def parse_grades(grades_name_dict: STR_DICT, header: int = 13, footer: int = 4) -> DF_DICT:
        df_dict = {}
        for name, file in grades_name_dict.items():
            dataframe = pd.read_fwf(open(file, 'r'), skiprows=header, skipfooter=footer)
            dataframe['ID'] = dataframe['3LATOM'].str.extract(pat=r"([A-Z]{3}\d+:[A-Z])")
            dataframe['RES'] = dataframe['ID'].str.extract(pat=r"([A-Z]{3})")
            dataframe['ID'] = dataframe['ID'].str.extract(pat=r"(\d+:[A-Z])")
            dataframe['ID'] = ':' + dataframe['ID'].str.replace(":", ".")
            df = dataframe[['RES', 'ID', 'SCORE']]
            all_graded = df.shape[0]
            df = df.dropna()
            on_structure = df.shape[0]
            # print(df.head(10))
            df_dict[name] = df
            print(f"{name}: {all_graded} graded residues, {on_structure} show up on structure")
            # print(df)
        return df_dict

    def merge_df_dict(asc_df_dict: DF_DICT, grades_df_dict: DF_DICT) -> pd.DataFrame:  # TODO: change for consurf
        hwy, hwx = [k for k in asc_df_dict.keys()]
        hwy_two, hwx_two = [k for k in grades_df_dict.keys()]
        if hwy != hwy_two or hwx != hwx_two:
            print(f"KEYS DO NOT MATCH UP: {hwy}:{hwy_two}, {hwx}:{hwx_two}")
            print(f"\nKeys:\n\tHWX: '{hwx}'\n\tHWY: '{hwy}'")
            sys.exit()
        print(f"\nKeys:\n\tHWX: '{hwx}'\n\tHWY: '{hwy}'")
        grades_asc_df_dict = {}
        for name in asc_df_dict.keys():
            df = asc_df_dict[name].merge(grades_df_dict[name], on=['ID', 'RES'])
            # print(df)
            grades_asc_df_dict[name] = df
        merge_df = grades_asc_df_dict[hwx].merge(grades_asc_df_dict[hwy], on='POS', suffixes=(f'_{hwx}', f'_{hwy}'))
        merge_df['DIFF'] = pd.Series(merge_df[f'SCORE_{hwx}'] - merge_df[f'SCORE_{hwy}']).abs()
        return merge_df

    def main_df_generation(asc_file: str, split_line: int, grade_names_dict: STR_DICT) -> pd.DataFrame:  # TODO: change for consurf
        asc_dataframe_dict = parse_three_letter_asc(asc_file, split_line)
        grades_dataframe_dict = parse_grades(grade_names_dict)
        merge_df = merge_df_dict(asc_dataframe_dict, grades_dataframe_dict)
        return merge_df

    df = main_df_generation(asc_file, asc_split_line, grades_dict)
    print(df.columns)
    fig, ax = plt.subplots(facecolor='w')
    x, y = '2HWX', '2HWY'
    print('\n\n\n\n', df)
    print(df['DIFF'].max())
    ax.scatter(df[f'SCORE_{x}'], df[f'SCORE_{y}'],
               c=df['DIFF'],
               cmap='Spectral',
               alpha=0.75,
               )
    # TODO: Add annotations! Maybe only add them to the ones above some threshold
    # line = np.linspace(-5, 5, 11)
    # plt.plot(line, line, 'k-', alpha=0.25)
    ax.set_xlabel(f'Consurf score of residues in {x}')
    ax.set_ylabel(f'Consurf score of residues in {y}')
    plt.show()
    return df


if __name__ == '__main__':
    # example()
    asc = r"Chimera_Files/200522_Match2HWXand2HWY.asc"
    grades_file_dict = {"2HWY": r"./SMG5_Consurf_Outputs/consurf.grades",
                        "2HWX": r"./SMG6_Consurf_Outputs/consurf.grades",
                        }
    # df = my_denovo_hydro(asc, 118, 'kd')
    my_denovo_consurf(grades_file_dict, asc, 118)
