-- Revert timeline_analysis_project_event_type_values

BEGIN;

DELETE FROM timeline.analysis_project_event_type WHERE id = 'config_added';
DELETE FROM timeline.analysis_project_event_type WHERE id = 'status_changed';
DELETE FROM timeline.analysis_project_event_type WHERE id = 'model_created';
DELETE FROM timeline.analysis_project_event_type WHERE id = 'cle_changed';
DELETE FROM timeline.analysis_project_event_type WHERE id = 'instrument_data_assigned';

COMMIT;
