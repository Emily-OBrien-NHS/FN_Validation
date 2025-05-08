import pandas as pd
import matplotlib.pyplot as plt

events_path = 'G:/PerfInfo/Performance Management/PIT Adhocs/2024-2025/HannahP 2425/UEC Adhocs/Frazer Nash UEC/Output/Run 6/ED Baseline with resources 01042024 28 days [1]/'
events = pd.read_csv(events_path + 'Events.csv')
#Filter to resource events and change to date time
events = events.loc[events['Event_Type'].isin(['Resource Assigned', 'Resource Released'])].copy()
events['DateTime'] = pd.to_datetime(events['DateTime'])

#Get each individual staff member as it's own row
events['Staff Member'] = events['Details'].str.split(', ')
events = events.explode('Staff Member')
events['Staff Type'] = (events['Staff Member'].str.split(' ').str[:-1]
                        .apply(lambda x: ' '.join(x)))
#Try to get the time taken for each task
events['Time'] = (events.sort_values(by=['DateTime', 'Event_Type',
                                         'Staff Member', 'Process', 'Details'])
                        .groupby(['Staff Member', 'Process', 'Node', 'Details'])
                                ['DateTime'].diff().dt.seconds / 60).to_list()
#Get dummies to get the number of times an event happens
dummies = pd.get_dummies(events['Process'])
events = events.join(dummies)
events['Date'] = pd.to_datetime(events['DateTime']).dt.date.astype(str)
events['Hour'] = pd.to_datetime(events['DateTime']).dt.hour.astype(int)
events['Date Hour'] = events['Date'] + ' ' + events['Hour'].astype(str)

#Get average number of tasks done hourly by each staff type
average_hourly = events.groupby(['Staff Member', 'Date Hour', 'Hour'], as_index=False)[dummies.columns].sum()
average_hourly['Staff'] = average_hourly['Staff Member'].str.split(' ').str[:-1].apply(lambda x: ' '.join(x))
average_hourly = average_hourly.groupby(['Staff', 'Hour'], as_index=False)[dummies.columns].mean().sort_values(by=['Staff', 'Hour'])

#Average time staff members spend on tasks
average_time = events.groupby(['Staff Member', 'Process', 'Hour'],
                              as_index=False)['Time'].mean()
average_time['Staff'] = (average_time['Staff Member'].str.split(' ').str[:-1]
                        .apply(lambda x: ' '.join(x)))
average_staff_time = average_time.groupby(['Staff', 'Process', 'Hour'],
                                          as_index=False)['Time'].mean().sort_values(by=['Staff', 'Process'])

#Write to excel
writer = pd.ExcelWriter('/'.join(events_path.split('/')[:-2])+'/Staff Usage.xlsx',
                        engine='xlsxwriter')   
average_hourly.to_excel(writer, sheet_name='Average Hourly', index=False)
average_staff_time.to_excel(writer, sheet_name='Average Time', index=False)
events.to_excel(writer, sheet_name='Events', index=False)
writer.close()



#Majors HCAs
Majors_HCAS_average_hourly = (average_hourly.loc[
                             average_hourly['Staff'].str.contains("HCA (Majors)",
                                                                  regex=False)]
                             .set_index(['Staff', 'Hour'])
                             .sort_values(by='Hour'))
Majors_HCAS_average_hourly[[col for col in Majors_HCAS_average_hourly.columns
                            if Majors_HCAS_average_hourly.sum()[col] > 0]]

Majors_HCAS_average_time = average_staff_time.loc[average_staff_time['Staff'].str.contains("HCA (Majors)",
                                                                  regex=False)].pivot(index=['Staff', 'Hour'], columns='Process', values='Time').plot()