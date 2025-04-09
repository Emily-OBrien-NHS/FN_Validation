import pandas as pd

events = pd.read_csv('G:/PerfInfo/Performance Management/PIT Adhocs/2024-2025/HannahP 2425/UEC Adhocs/Frazer Nash UEC/Output/Run 4/ED Baseline with resources 01042024 28 days [1]/Events.csv')
events = events.loc[events['Event_Type'] == 'Resource Assigned'].copy()
events['Staff Member'] = events['Details'].str.split(', ')
events = events.explode('Staff Member')
events = events.join(pd.get_dummies(events['Process']))

#Majors HCAs
Majors_HCAs = events.loc[events['Staff Member'].str.contains("HCA (Majors)", regex=False),
['DateTime', 'Staff Member', 'Admitted (Majors)', 'Bloods/ECG (Majors)',
'Misc Assessment (Majors)', 'Nurse Assessment (Majors)', 'Obs 30 min (Majors)',
'Obs 60 min (Majors)', 'Treatment (Majors)', 'Triaged (Majors)']].copy()
Majors_HCAs['Hour'] = pd.to_datetime(Majors_HCAs['DateTime']).dt.date.astype(str) + ' ' +  pd.to_datetime(Majors_HCAs['DateTime']).dt.hour.astype(str)
Majors_HCAs.groupby(['Staff Member', 'Hour'], as_index=False).sum().groupby('Hour').mean(numeric_only=True).plot()
