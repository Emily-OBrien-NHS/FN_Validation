import data_cleaning_and_transformation as cleaning
import pathway_definitions as pathways
import process_durations as durations
import validation
from sqlalchemy import create_engine
import config
import pandas as pd
import os

if __name__ == "__main__":
    #Set up folder to store this scenario and other subfoleders on the shared
    # drive from filepath in config file
    if not os.path.exists(config.output_path):
      os.mkdir(config.output_path)
    os.chdir(config.output_path)
    #If no plots folder in this path, add it in
    other_dirs = ["Additional Outputs",
                  "Additional Outputs/Duration Distributions",
                  "Additional Outputs/Flow Diagrams"]
    for dir in other_dirs:
        if not os.path.exists(dir):
            os.mkdir(dir)

################################################################################
#-------------------------Read in and clense raw data--------------------------#
################################################################################
    start_date = '01-APR-2024 00:00:00'
    end_date = '30-APR-2024 23:59:59'
    # ---------------------- Database Engines
    realtime_engine = create_engine('mssql+pyodbc://@dwrealtime/RealTimeReporting?'\
                           'trusted_connection=yes&driver=ODBC+Driver+17'\
                               '+for+SQL+Server')
    cl3_engine = create_engine('mssql+pyodbc://@cl3-data/DataWarehouse?'\
                           'trusted_connection=yes&driver=ODBC+Driver+17'\
                               '+for+SQL+Server')
    # ---------------------- Events data
    events_query = f"""SET NOCOUNT ON
                    -------CONNECT TO DWREALTIME

                    -------First get attendance-level data
                    select AttendanceID = ncattendanceId
                            ,ArrivalMode = case when AmbulanceArrivalDate is NULL then 'Walk-In' else 'Ambulance' end
                            ,Injury = IsInjury
                            ,TriageCategory
                            ,ArrivalDateTime
                            ,DischargeDateTime
                    into #att
                    from [cl3-data].DataWarehouse.ed.vw_EDAttendance
                    where dischargedatetime between '{start_date}' and '{end_date}'

                    ----Use location table to get all locations for these attendances
                    select NCAttendanceId, LocationDescription, LocationSubType, LocationOrder, StartDateTime, EndDateTime
                    into #locs
                    from [cl3-data].DataWarehouse.ed.vw_EDAttendanceLocationHistory loc 
                    inner join #att att on att.AttendanceID = loc.NCAttendanceId

                    --------------------------------------------------------------------------------------------------------------------
                    ---------------ADMISSIONS QUERY-------------------------------------------------------------------------------------
                    --------------------------------------------------------------------------------------------------------------------

                    ----------Find which patients go on to be admitted from the attendances
                    select nerve.NCAttendanceID, case when admitprvsprefno is not NULL and ActualDischargeDestinationWardCode in ('rk950aau','rk950aau01', 'rk950afu') then 'Admitted - SDEC'
                                                when admitprvsprefno is not NULL and ActualDischargeDestinationWardCode in ('rk950mau','rk950amw') then 'Admitted - MAU'
                                                when admitprvsprefno is not NULL and ActualDischargeDestinationWardCode like 'rk950%' then 'Admitted - Other Derriford Ward'
                                                else 'Discharged' end as Adm,
                                                ActualDischargeDestinationWardCode
                    into #adm
                    from [cl3-data].DataWarehouse.ed.vw_EDAttendance nerve
                    ---inner join attendances to get only required records
                    inner join #att att on att.AttendanceID = nerve.NCAttendanceId

                    ----Uncomment this bit to get table of where patients arrived to
                    --select att.*, 
                    --case when LocationSubType = 'Waiting Area' then 'Ambulatory Waiting Area'
                    --	when LocationSubType = 'Amb Cubicles' then 'Ambulatory Cubicles'
                    --	when LocationSubType = 'Minors Paeds' then 'Paediatrics'
                    --	when LocationSubType = 'Ambulance Bay/HALO' then 'Ambulance'
                    --	when LocationSubType = 'Corridor' then 'Majors Corridor'
                    --	else LocationSubType end as InitialLocation
                    --from #att att
                    --left join #locs locs on locs.NCAttendanceId = att.AttendanceID
                    --				and LocationOrder = 1 -----Only get the initial location on arrival

                    ---Before getting the events table, need to pull through
                    ----first clinician contact - agreed with Nanette to choose the first of either seen by, clerking or senior review 
                    select VisitId, min([timestamp]) as timestamp--, ArrivalDateTime, DischargeDateTime
                    into #seenbytime --drop table #events
                    FROM [NerveCentreFeed].[Note].[ClinicalNoteHistory] note
                    inner join #att att on att.AttendanceID = note.VisitId ---only get attendances within timeframe
                    where strikeoutid is NULL ---don't pick up cancelled/overwritten values
                    and NoteKey in ('ED Seen By','ED Senior Reviewed',
                        'ED Clerking Actual Date/Time')
                    group by VisitId


                    ---Also need to get the first bloods/ecg for each attendance
                    select NCattendanceID, min(DateRaised) as timestamp--, ArrivalDateTime, DischargeDateTime
                    into #bloodstime --drop table #events
                    FROM [NerveCentreFeed].[ED].[vw_EDTask] task
                    inner join #att att on att.AttendanceID = task.NCattendanceID ---only get attendances within timeframe
                    where [CloseReason] = 'Completed'
                    and TaskCategory in ('ED Task 01 - Bloods','ED Task 04 - ECG')
                    group by NCattendanceID

                    --------------------------------------------------------------------------------------------------------------------
                    ---------------MAIN EVENTS QUERY------------------------------------------------------------------------------------
                    --------------------------------------------------------------------------------------------------------------------

                    -----Get event table
                    ---First get a temp table from the note table
                    select VisitId,NoteKey, NoteValue, [timestamp],AddedBy, ArrivalDateTime, DischargeDateTime
                    into #events --drop table #events
                    FROM [NerveCentreFeed].[Note].[ClinicalNoteHistory] note
                    inner join #att att on att.AttendanceID = note.VisitId ---only get attendances within timeframe
                    where strikeoutid is NULL ---don't pick up cancelled/overwritten values
                    and NoteKey in ('ED Arrival Transport Mode','ED Ambulance Arrival Date/Time'--,'ED Seen By',
                        ,'ED Nurse completing Triage','Manchester Triage Score', 'ED Senior Reviewed'
                        --,'ED Clerking Actual Date/Time'
                        ,'ED Specialty Reviewed Dt','DTA Actual Date/Time','ED Departure Ready Date/Time','ED Discharge Clinician')

                    ---union nursing assessments data
                    union all
                    ---Get nursing assessmnets and group them up
                    select VisitId
                        ,case when NoteKey in ('Sepsis Screening', 'Falls Risk','Mental Health Risk','Patient has capacity','ED Pain') then 'Misc Assessments' else NoteKey end as NoteKey
                        ,case when NoteKey in ('Sepsis Screening', 'Falls Risk','Mental Health Risk','Patient has capacity','ED Pain') then 'Misc Assessments' else NoteKey end as NoteValue
                        , min([timestamp]) as [timestamp]--
                        ,AddedBy  = NULL
                        , att.ArrivalDateTime, att.DischargeDateTime 
                    --into #nurs_assess
                    from  [NerveCentreFeed].[Note].[ClinicalNoteHistory] note
                    inner join #att att on att.AttendanceID = note.VisitId ---only get attendances within timeframe
                    where strikeoutid is NULL ---don't pick up cancelled/overwritten values
                    and notekey in ('ED Nursing Assessment By' --Nursing assessment
                                    ,'Sepsis Screening', 'Falls Risk','Mental Health Risk','Patient has capacity','ED Pain') --Misc Assessments
                                    --,'ED Secondary Nursing Assessment')--Nursing continuations
                    and NoteValue <> 'Not Screened for Sepsis' --This is autofilled and gets updated later
                    group by
                    VisitId
                    , att.ArrivalDateTime, att.DischargeDateTime 
                    ,case when NoteKey in ('Sepsis Screening', 'Falls Risk','Mental Health Risk','Patient has capacity','ED Pain') then 'Misc Assessments' else NoteKey end

                    ----also union new version of seen by - only using first value
                    union all
                    select sbt.VisitId, NoteKey = 'ED Seen By',NoteValue = 'ED Seen By', sbt.[timestamp],  addedby,ArrivalDateTime, DischargeDateTime
                    from #seenbytime sbt
                    left join [NerveCentreFeed].[Note].[ClinicalNoteHistory] note
                            on note.VisitId = sbt.VisitId
                            and note.NoteKey in ('ED Seen By','ED Senior Reviewed',
                        'ED Clerking Actual Date/Time')
                        and sbt.timestamp = note.timestamp
                    inner join #att att on att.AttendanceID = note.VisitId ---only get attendances within timeframe

                    --also add in treatment data
                    union all
                    ----------------Treatments from tasks
                    select NCattendanceID as VisitID
                            ,'Treatment'
                            ,TaskCategory
                            ,[CompletedDateTime]
                            --,TreatmentNumber = 'Treatment' + cast(ROW_NUMBER () over(partition by NCattendanceID order by [CompletedDateTime] asc) as varchar
                            ,AddedBY = [ClosedByUsername]
                            ,att.ArrivalDateTime, att.DischargeDateTime 
                    FROM [NerveCentreFeed].[ED].[vw_EDTask] task
                    inner join #att att on att.AttendanceID = task.ncattendanceID ---only get attendances within timeframe
                    where [CloseReason] = 'Completed'
                    ----Only use ED tasks that are treatments, as discussed with Nanette
                    and TaskCategory in ('ED Task 02 - Cannula'
                                            ,'ED Task 07 - ePMA Prescription'
                                            ,'ED Task 08 - ePMA Critical Meds'
                                            ,'ED Task 10 - Plaster of Paris'
                                            ,'ED Task 11 - Minor Injury Treatment'
                                            ,'ED Task 15 - Urinary Catheter'
                                            ,'ED Task Ametop'--space at the end?
                                            ,'ED Task Blood Transfusion'
                                            ,'ED Task Cervical Collar Placement'
                                            ,'ED Task Fascia Iliaca Block'
                                            ,'ED Task IV Fluids'
                                            ,'ED Task Manipulation required'
                                            ,'ED Task Medication Required'
                                            ,'ED Task NG Tube Placement'
                                            ,'ED Task Removal of Plaster of Paris'
                                            ,'ED Task Removal of Ring'
                                            ,'ED Task Splint/Strap/Sling'
                                            ,'ED Task Split/Strap/Sling'
                                            ,'ED Task Suture'
                                            ,'ED Task TTAs'
                                            ,'ED Task TWOC'
                                            ,'ED Task Wound Care')


                    ----also add in bloods/ECG
                    union all


                    select task.NCattendanceID as VisitID
                            ,'Bloods/ECG'
                            ,TaskCategory
                            ,task.timestamp
                            ,AddedBy = [ClosedByUsername]
                            ,att.ArrivalDateTime, att.DischargeDateTime 
                            FROM #bloodstime task
                    inner join #att att on att.AttendanceID = task.ncattendanceID ---only get attendances within timeframe
                    left join [NerveCentreFeed].[ED].[vw_EDTask] tsk on tsk.NCAttendanceId = task.NCAttendanceId
                                and tsk.DateRaised = task.timestamp
                                and tsk.TaskCategory in ('ED Task 01 - Bloods','ED Task 04 - ECG')



                    --------------------------------------------------------------------------------------------------------------------
                    ---------------FORMATTED EVENTS OUTPUT------------------------------------------------------------------------------
                    --------------------------------------------------------------------------------------------------------------------

                    select VisitId, 
                    ---Next, get the event name
                    case when NoteKey = 'ED Arrival Transport Mode' then 'Booked In'
                        when NoteKey = 'ED Ambulance Arrival Date/Time' then 'Ambulance Arrival'
                        when NoteKey = 'ED Seen By' then 'Seen By Clinician/Treated'
                        when NoteKey = 'Treatment' then 'Treatment'
                        when NoteKey = 'Bloods/ECG' then 'Bloods/ECG'
                        when NoteKey = 'ED Nursing Assessment By' then 'Nurse Assessment'
                        when NoteKey = 'Misc Assessments' then 'Misc Assessment'
                        when NoteKey in ('ED Nurse completing Triage','Manchester Triage Score') then 'Triaged'
                        when NoteKey = 'ED Senior Reviewed' then 'Senior Reviewed'
                        when NoteKey = 'ED Clerking Actual Date/Time' then 'Clerked'
                        when NoteKey = 'ED Specialty Reviewed Dt' then 'Specialty Reviewed'
                        when NoteKey = 'DTA Actual Date/Time' then 'Decision to Admit'
                        when NoteKey = 'ED Departure Ready Date/Time' then 'Clinically Ready to Proceed'
                        when NoteKey = 'ED Discharge Clinician' then adm.Adm ---For the discharge event, use formatted admission details from temp table
                        end as EventName
                    ----Then get the timestamp
                    ,case when NoteKey in ('ED Ambulance Arrival Date/Time','ED Clerking Actual Date/Time','ED Specialty Reviewed Dt','ED Departure Ready Date/Time')
                            then [NerveCentreFeed].Util.ConvertNCDateNumericToLocalDateTime(NoteValue)
                            when NoteKey = 'ED Arrival Transport Mode' then ArrivalDateTime
                            when NoteKey = 'ED Discharge Clinician' then DischargeDateTime
                        else [timestamp] end as EventTime
                    ----Get staff member
                    ,case when NoteKey in ('ED Nurse completing Triage'--,'ED Seen By'
                                        ,'ED Discharge Clinician')
                                then NoteValue
                        else AddedBy end as EventStaffMember
                    ----Replace staff names with numbers
                    ,dense_rank() over  (order by case when NoteKey in ('ED Nurse completing Triage','ED Seen By','ED Nursing Assessment By','ED Discharge Clinician')
                                then NoteValue
                        else AddedBy end) as EventStaffId
                    ,case when LocationSubType = 'Waiting Area' then 'Ambulatory Waiting Area'
                        when LocationSubType = 'Amb Cubicles' then 'Ambulatory Cubicles'
                        when LocationSubType = 'Minors Paeds' then 'Paediatrics'
                        when LocationSubType = 'Ambulance Bay/HALO' then 'Ambulance'
                        when LocationSubType = 'Corridor' then 'Majors Corridor'
                        else LocationSubType end as EventLocation	  
                    --,NoteKey, NoteValue, [timestamp],AddedBy
                    FROM #events note--[NerveCentreFeed].[Note].[ClinicalNoteHistory] note
                    ---join admission status
                    left join #adm adm on adm.ncattendanceID = note.VisitId
                    --inner join #att att on att.AttendanceID = note.VisitId ---only get attendances within timeframe
                    left join #locs locs on locs.NCAttendanceId = note.VisitId ---get location at event time
                            and case when NoteKey in ('ED Ambulance Arrival Date/Time','ED Clerking Actual Date/Time','ED Specialty Reviewed Dt','ED Departure Ready Date/Time')
                            then [NerveCentreFeed].Util.ConvertNCDateNumericToLocalDateTime(NoteValue)
                            when NoteKey = 'ED Arrival Transport Mode' then note.ArrivalDateTime
                            when NoteKey = 'ED Discharge Clinician' then note.DischargeDateTime
                        else [timestamp] end between locs.StartDateTime and locs.EndDateTime
                    --where VisitId = '1395736'
                    --and strikeoutid is NULL ---don't pick up cancelled/overwritten values
                    -- and NoteKey in ('ED Arrival Transport Mode','ED Ambulance Arrival Date/Time','ED Seen By',
                        --'ED Nursing Assessment By','ED Nurse completing Triage','Manchester Triage Score', 'ED Senior Reviewed',
                        --'ED Clerking Actual Date/Time','ED Specialty Reviewed Dt','DTA Actual Date/Time','ED Departure Ready Date/Time','ED Discharge Clinician')
                    order by VisitId




                    ------------------Additional asks - Obs
                    --select NCAttendanceId, ChartType, ChartDateTime
                    --	into #obs --select * from #obs
                    --	 FROM [cl3-data].[DataWarehouse].[ED].[vw_EDAttendanceObservationChartTotal]
                    --	 where ChartDateTime between '25-MAR-2023 00:00:00' and '31-MAR-2024 23:59:59'

                    ----Get obs that correspond to attendances provided
                    --select obs.*
                    --from #obs obs
                    --inner join #att att on att.AttendanceID = obs.NCAttendanceId
                    --order by obs.NCAttendanceID



                    ----------------------------------------------------------------------------------------------------------------------
                    -----------------ADMISSIONS QUERY-------------------------------------------------------------------------------------
                    ----------------------------------------------------------------------------------------------------------------------

                    ------------Find which patients go on to be admitted from the attendances
                    --select nerve.NCAttendanceID, case when admitprvsprefno is not NULL and ActualDischargeDestinationWardCode in ('rk950aau','rk950aau01', 'rk950afu') then 'Admitted - SDEC'
                    --							when admitprvsprefno is not NULL and ActualDischargeDestinationWardCode in ('rk950mau','rk950amw') then 'Admitted - MAU'
                    --							when admitprvsprefno is not NULL and ActualDischargeDestinationWardCode like 'rk950%' then 'Admitted - Other Derriford Ward'
                    --							else 'Non-Admitted' end as Adm,
                    --							ActualDischargeDestinationWardCode
                    --into #adm
                    --from [cl3-data].DataWarehouse.ed.vw_EDAttendance nerve
                    -----inner join attendances to get only required records
                    --inner join #att att on att.AttendanceID = nerve.NCAttendanceId
                    """
    events_raw = pd.read_sql(events_query, realtime_engine)
    # ---------------------- Imaging data
    imaging_query = f"""SET NOCOUNT ON
    declare @startdttm as datetime
 declare @enddttm as datetime

 set @startdttm = '{start_date}'
 set @enddttm = '{end_date}'

--Get CRiS Reports
  select EventKey, HospitalNumber,
  case when max(SubModality) = 'Radiology' then max(isnull(EventDateTime,'1900-01-01')) else max(isnull(VerifiedDateTime,'1900-01-01')) end as ReportedDateTime,--If XR, use event time, otherwise use report verified
  max(SubModality) as SubModality
  into  #cris_rep
  from [DataWarehouse].[Imaging].[vw_Report]
  where SubModality in ('CT','MRI','Radiology','Ultrasound') --Only get imaging scans
  and PatientTypeCode in ('C','J') --only ED patients
  and CreationDateTime between @startdttm and @enddttm
  and SiteCode = 'RK950' --Don't include MIU/UTC/Community scans
  group by EventKey, HospitalNumber


--GET CRIS IMAGING DATA
  select cris.EventKey, cris.HospitalNumber, min(isnull(CreationDateTime,'2099-01-01')) as OrderEnteredDateTime,
  case when (max(isnull(cris_rep.ReportedDateTime,'1900-01-01')) is NULL or max(isnull(cris_rep.ReportedDateTime,'1900-01-01')) = '1900-01-01') and max(cris.SubModality) = 'Radiology' then max(isnull(cris.EventDateTime,'1900-01-01')) else max(isnull(cris_rep.ReportedDateTime,'1900-01-01')) end as ResultsAvailableDateTime,
  max(cris.SubModality) as TestName,
  ItemMasterCategory = 'Imaging'
  ,UrgencyCode = cris.UrgencyCode, Urgency = cris.Urgency
  into #cris --drop table #cris
  FROM [DataWarehouse].[Imaging].[vw_Activity] cris
  --join reports
  left join #cris_rep cris_rep on cris_rep.EventKey = cris.EventKey
  where cris.SubModality in ('CT','MRI','Radiology','Ultrasound') --Only get imaging scans
  and cris.StatusCode = 'ATP' --remove those who DNA/not performed
  and cris.PatientTypeCode in ('C','J') --only ED patients
  and cris.CreationDateTime between @startdttm and @enddttm
  and SiteCode = 'RK950' --Don't include MIU/UTC/Community scans
  group by cris.EventKey, cris.HospitalNumber,cris.UrgencyCode,cris.Urgency


--Get temp ED data
select nerve.NCAttendanceId
              ,nerve.HospitalNumber
              ,nerve.ArrivalDateTime
			  ,nerve.DischargeDateTime
              ,nerve.IsNewAttendance
into #nerve
from   DataWarehouse.ed.vw_EDAttendance nerve
where   nerve.dischargedatetime between @startdttm and @enddttm
and nerve.DischargeDateTime is not NULL



---JOIN TEST REQUESTS TO ED ATTENDS 
select nerve.NCAttendanceId AS VisitId         
,nerve.ArrivalDateTime               
,res.OrderEnteredDateTime AS EventTime
,res.ResultsAvailableDateTime
,res.TestName AS EventName
from   #nerve nerve 
inner   join #cris res on nerve.HospitalNumber=res.HospitalNumber   --inner join so we don't return patients with no tests 
and           res.OrderEnteredDateTime between nerve.ArrivalDateTime and nerve.DischargeDateTime 
where  nerve.IsNewAttendance = 'y'  
and ResultsAvailableDateTime > '30-MAR-2024' ---exclude records with missing data
group by nerve.NCAttendanceId               
,nerve.ArrivalDateTime               
,res.OrderEnteredDateTime  
,res.ResultsAvailableDateTime
,res.TestName


    """
    imaging_raw = pd.read_sql(imaging_query, cl3_engine)
    # ---------------------- Close Connections
    realtime_engine.dispose()
    cl3_engine.dispose()
    
    # ---------------------- Events data
    #Add imaging data to events file
    imaging_events = imaging_raw[["VisitId", "EventName"]].copy()
    first_event_time = imaging_raw.groupby("VisitId",
                                           as_index=False)['EventTime'].min()
    imaging_events = imaging_events.merge(first_event_time, on='VisitId',
                                          how='left')
    events_raw = pd.concat([events_raw, imaging_events])
    #drop duplicates and fix anomaly times
    events_quality = cleaning.drop_duplicates_and_anomaly_times_events_data(
                     events_raw, config.repeat_time_threshold,
                     config.remove_duplicate_staffid,
                     config.remove_duplicate_location)

    # ---------------------- Combine and Clense data
    events_quality, treat_repeat = cleaning.main_cleanse_and_transform_data(
                                   events_quality,
                                   config.event_names_to_exclude_for_repetition,
                                   config.excluded_event_names,
                                   config.locations_to_drop,
                                   config.natural_order_for_processes,
                                   config.include_spawn_end_events,
                                   config.locations_pathway_map,
                                   config.keep_last_location,
                                   config.admitted_map)

################################################################################
#-----------------------------Pathway Definitions------------------------------#
################################################################################

    # ---------------------- Calculate the number of patients making each
    #                        transition.
    transitions = pathways.add_reset_transitions(events_quality)

    #filters out any transitions below % threshold specified in config file if
    #bools are set toFale, True.  If both False, no filtering occurs, and if
    #True, False then any patient who has a transition below Threshold%
    # is removed.
    pathways.main_generate_dfg_and_pathway_definitions(transitions,
                 config.output_path, config.include_spawn_end_events,
                 config.obs_splits, config.export_event_log_csv,
                 config.export_log_to_csv_after_using_log_converter,
                 config.pathways_wait_in_place, config.transition_threshold, 
                 False, True, "EventName", "Pathway")


################################################################################
#------------------------------Process Durations-------------------------------#
################################################################################

    # ------------------------ Calculate the duration of each event in the data
    # ------------------------ and create histograms and process durations
    event_diffs = durations.add_difference_in_minutes_to_durations(
                  events_quality, config.where_duration_should_be_0)
    event_diffs = durations.add_imaging_timings(event_diffs, imaging_raw)

    #Max threshold 14 hours, and config.quantile_threshold percentile"
    durations.main_generate_histogram_and_process_durations(
              event_diffs, treat_repeat, config.output_path,
              [durations.within_threshold_diff(840),
               durations.within_diff_quantile(config.quantile_threshold)])
    
    # ------------------------Other filtering options, replace list of functions
    # -- durations for different filterings/scenarios
    # analysis_name = "Max threshold 2 hours and including 100 percentile"
    # [durations.within_threshold_diff(120)]
    # -- analysis_name = "Max threshold 2 hours, 100 perc, between 8am and 10pm"
    # [durations.within_threshold_diff(120), durations.only_daytime_events]
    # -- analysis_name = f"Max threshold 14 hours,
    # --and {int(config.quantile_threshold*100)} percentile"
    # [durations.within_threshold_diff(840),
    #  durations.within_diff_quantile(config.quantile_threshold)]
    # -- analysis_name = f"Max threshold 14 hours,
    # --{int(config.quantile_threshold*100)} perc, between 8am and 10pm"
    # [durations.within_threshold_diff(840),
    #  durations.within_diff_quantile(config.quantile_threshold),
    #  durations.only_daytime_events]


################################################################################
#------------------------Add other files for scenario--------------------------#
################################################################################
#Read in data from existing input file and save into this location to create a
#full scenario
for file in config.additional_filenames:
    df = pd.read_csv(config.other_input_filepath + '/' + file + '.csv')
    df.to_csv(config.output_path + '/' + file + '.csv', index=False)

################################################################################
#---------------------------Add Summary Text File------------------------------#
################################################################################
#Get all the variables from the config file
config_variables = {name: value 
                    for name, value in vars(config).items() 
                    if not name.startswith('__')}
config_variables = {'start_date':start_date, 'end_date':end_date} + config_variables
#Create string of variables and their values at run time
output_str = ''
for name, value in config_variables.items():
    variable_str = str(name) + ' = ' + str(value)
    output_str = validation.print_and_add_str(output_str, variable_str)
#Save as txt file
config_output_path = config.output_path + '/Additional Outputs/Config Values.txt'
with open(config_output_path, 'w', encoding='utf-8') as f:
        f.write(output_str)

################################################################################
#-------------------------------Run Validation---------------------------------#
################################################################################
validation.validation_process(config.output_path)