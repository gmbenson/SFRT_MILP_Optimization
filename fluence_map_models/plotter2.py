import matplotlib

matplotlib.use('TkAgg')  # Fix backend issues in PyCharm

import pandas as pd
import matplotlib.pyplot as plt
import seaborn as sns
from io import StringIO
from matplotlib.ticker import MaxNLocator



# Updated dataset with patient 21
data = """
patient_id	NCI CT	ptv volume	Candidate Spheres	PTV voxels	total voxels	Num Nodes	Gap	GTV Dev	PTV Dev	Esophagus Dev	Heart Dev	Left Lung Dev	Right Lung Dev	spinal cord dev	skin dev	body dev	Total Dev	Average Dev	Time
4	5	312.6	273	5298	159110	10	0.6326	14.85	16.32	13.99	0.79	11.98	5.08	18.97	10.36	8.98	1433460.902	9.009244562	40hr
4	5	312.6	273	5298	159110	9	0.6753	14.81	16.3	14.04	0.79	11.98	5.07	18.98	10.36	8.98	1433854.578	9.0117188	40hr
4	5	312.6	273	5298	159110	8	0.7159	14.87	16.36	14.05	0.79	11.98	5.06	18.98	10.37	8.99	1434530.587	9.01596749	40hr
14	22	357.2	371	5930	166986	10	0.7713	7.87	7.19	8.95	11.21	6.16	8.17	23.03	7.66	7.75	1293462.175	7.745931843	40hr
14	22	357.2	371	5930	166986	9	0.8686	7.43	7.01	8.92	11.2	6.14	8.17	23.03	7.66	7.74	1292412.065	7.739643232	40hr
14	22	357.2	371	5930	166986	8	0.8198	7.34	6.9	8.9	11.22	6.1	8.15	23.03	7.66	7.74	1291573.307	7.73462031	40hr
15	23	1276.9	2213	19769	212366	42	2.0098	4.57	6.32	7.89	9.29	13.7	9.11	35.33	12.29	10.68	2268576.719	10.68239134	40hr
15	23	1276.9	2213	19769	212366	41	2.0096	4.56	6.26	8.01	9.26	13.72	9.08	35.33	12.3	10.68	2268352.148	10.68133386	40hr
15	23	1276.9	2213	19769	212366	39	1.9832	4.43	6.21	7.94	9.27	13.69	9.07	35.35	12.32	10.67	2267396.227	10.67683258	40hr
17	25	372.6	360	6177	125451	11	0	5.48	7.43	8.49		7.62	5.68	16.28	5.93	5.29	665618.7669	5.305806784	40hr
17	25	372.6	360	6177	125451	10	1.7454	5.1	7.3	8.63		7.62	5.69	16.21	5.95	5.29	665729.8462	5.306692224	40hr
17	25	372.6	360	6177	125451	9	1.6543	4.58	7.03	8.64		7.63	5.68	16.18	5.95	5.27	664085.2592	5.293582827	40hr
21	33	662.6	538	10818	185072	17	1.253	8.56	9.8	10.25	5.82	8.56	7.49	15.59	10.13	8.07	1496651.326	8.08685985	40hr
21	33	662.6	538	10818	185072	16	1.2793	8.44	9.73	10.26	5.82	8.56	7.48	15.6	10.14	8.07	1496014	8.083416183	40hr
21	33	662.6	538	10818	185072	15	1.25	8.36	9.67	10.26	5.82	8.56	7.48	15.6	10.13	8.06	1495439.964	8.080314493	40hr
22	34	348.3	334	5739	136503	10	0.7463	7.83	11.24	6.29	1.32	6.82	8.62	30.87	6.29	6.96	950640.8505	6.964248775	40hr
22	34	348.3	334	5739	136503	9	1.1578	7.57	10.92	6.44	1.32	6.81	8.61	30.8	6.31	6.96	950174.326	6.960831088	40hr
22	34	348.3	334	5739	136503	8	1.0512	6.73	10.77	6.3	1.32	6.75	8.63	30.87	6.3	6.95	949034.3001	6.952479433	40hr
24	40	733.2	908	11717	191620	21	1.8096	5.57	6.04	13.9	9.74	12.31	5.67	24.88	8.89	7.82	1498624.483	7.820814542	40hr
24	40	733.2	908	11717	191620	20	1.8021	5.4	5.9	13.94	9.72	12.31	5.66	24.92	8.9	7.82	1497972.282	7.81741093	40hr
24	40	733.2	908	11717	191620	19	1.7419	5.32	5.84	13.95	9.71	12.31	5.65	24.91	8.9	7.81	1497409.278	7.814472801	40hr
25	42	476.3	580	7739	187217	14	1.2778	8.09	7.81	5.52	12.14	11.76	6.45	24.95	7.78	7.48	1400727.697	7.481840308	40hr
25	42	476.3	580	7739	187217	13	 1.2517 	7.98	7.74	5.52	12.14	11.76	6.45	24.95	7.78	7.48	1400189.351	7.478964792	40hr
25	42	476.3	580	7739	187217	12	 1.2140 	7.86	7.66	5.51	12.14	11.76	6.45	24.95	7.78	7.47	1399658.45	7.476129035	40hr

"""

# Load data
df = pd.read_csv(StringIO(data), sep='\t')

# Strip whitespace from string columns
df = df.applymap(lambda x: x.strip() if isinstance(x, str) else x)

# Convert relevant columns to numeric
for col in ['Average Dev', 'PTV Dev', 'Total Dev', 'Num Nodes']:
    df[col] = pd.to_numeric(df[col], errors='coerce')

# Drop rows with missing values in key columns
df = df.dropna(subset=["Average Dev", "PTV Dev", "Total Dev", "Time"])



cols_to_convert = [
    'Num Nodes', 'Total Dev', 'PTV Dev', 'GTV Dev',
    'Esophagus Dev', 'Heart Dev', 'Left Lung Dev', 'Right Lung Dev',
    'spinal cord dev', 'skin dev', 'body dev'
]
df[cols_to_convert] = df[cols_to_convert].apply(pd.to_numeric, errors='coerce')

# Drop rows with missing Num Nodes (required for x-axis)
df = df.dropna(subset=['Num Nodes'])

# Define unique patient IDs
patients = df['patient_id'].unique()

for patient in patients:
    patient_df = df[df['patient_id'] == patient].copy()
    patient_df = patient_df.sort_values('Num Nodes')

    # Plot 1: Total Dev
    plt.figure(figsize=(8, 5))
    ax = sns.lineplot(data=patient_df, x='Num Nodes', y='Total Dev', marker='o')
    ax.xaxis.set_major_locator(MaxNLocator(integer=True))
    print(f'Patient {patient}: Total Dev vs Num Nodes')
    plt.xlabel('Number of Nodes (MIS Size)')
    plt.ylabel('Total Deviation')
    plt.tight_layout()
    plt.show()

    # Plot 2: Scaled OAR Deviations
    oar_cols = ['Esophagus Dev', 'Heart Dev', 'Left Lung Dev', 'Right Lung Dev',
                'spinal cord dev', 'skin dev', 'body dev']
    scaled_oars = patient_df[oar_cols].copy()
    valid_max = scaled_oars.max().replace(0, pd.NA).dropna()
    scaled_oars = scaled_oars[valid_max.index].divide(valid_max)

    plt.figure(figsize=(10, 6))
    ax = plt.gca()
    for col in scaled_oars.columns:
        if scaled_oars[col].notna().any():
            plt.plot(patient_df['Num Nodes'], scaled_oars[col], marker='o', label=col)

    ax.xaxis.set_major_locator(MaxNLocator(integer=True))
    print(f'Patient {patient}: OAR Deviations (Normalized by Max) vs Num Nodes')
    plt.xlabel('Number of Nodes (MIS Size)')
    plt.ylabel('Deviation (Normalized)')
    plt.legend(title="OAR")
    plt.tight_layout()
    plt.show()

    # Plot 3: PTV & GTV
    plt.figure(figsize=(8, 5))
    ax = plt.gca()
    plt.plot(patient_df['Num Nodes'], patient_df['PTV Dev'], marker='o', label='PTV Dev')
    plt.plot(patient_df['Num Nodes'], patient_df['GTV Dev'], marker='o', label='GTV Dev')
    ax.xaxis.set_major_locator(MaxNLocator(integer=True))
    print(f'Patient {patient}: PTV & GTV Deviations vs Num Nodes')
    plt.xlabel('Number of Nodes (MIS Size)')
    plt.ylabel('Deviation')
    plt.legend()
    plt.tight_layout()
    plt.show()