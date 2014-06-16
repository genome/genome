-- Verify model.build_direct_status_columns

BEGIN;

SELECT status,date_scheduled,date_completed,run_by,created_by,created_at,updated_at FROM model.build WHERE FALSE;

SELECT 1/count(*) FROM pg_class WHERE relkind = 'i' and relname = 'build_status_idx';
SELECT 1/count(*) FROM pg_class WHERE relkind = 'i' and relname = 'build_date_scheduled_idx';
SELECT 1/count(*) FROM pg_class WHERE relkind = 'i' and relname = 'build_date_completed_idx';
SELECT 1/count(*) FROM pg_class WHERE relkind = 'i' and relname = 'build_run_by_idx';

ROLLBACK;
