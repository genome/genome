-- Revert timeline_allocation_event_type_values

BEGIN;

select 1/count(*) from timeline.analysis_project_event_type where id = 'config_added';
select 1/count(*) from timeline.analysis_project_event_type where id = 'status_changed';
select 1/count(*) from timeline.analysis_project_event_type where id = 'model_created';
select 1/count(*) from timeline.analysis_project_event_type where id = 'cle_changed';
select 1/count(*) from timeline.analysis_project_event_type where id = 'instrument_data_assigned';

COMMIT;
