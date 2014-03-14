-- Verify analysis_project_event_logging

BEGIN;
  SELECT id from timeline.analysis_project_event_type WHERE FALSE;
  SELECT status,is_cle,run_as from timeline.analysis_project WHERE FALSE;

  SELECT 1/count(*) FROM pg_class WHERE relkind = 'i' and relname = 'timeline_anp_object_id_index';

  SELECT 1/has_table_privilege('genome', 'timeline.analysis_project', 'TRUNCATE')::int;
  SELECT 1/has_table_privilege('gms-user', 'timeline.analysis_project', 'SELECT')::int;

  SELECT 1/has_table_privilege('genome', 'timeline.analysis_project_event_type', 'TRUNCATE')::int;
  SELECT 1/has_table_privilege('gms-user', 'timeline.analysis_project_event_type', 'SELECT')::int;

ROLLBACK;
