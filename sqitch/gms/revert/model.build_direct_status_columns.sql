-- Revert model.build_direct_status_columns

BEGIN;

DROP INDEX model.build_status_idx;
DROP INDEX model.build_date_scheduled_idx;
DROP INDEX model.build_date_completed_idx;
DROP INDEX model.build_run_by_idx;

ALTER TABLE model.build DROP COLUMN status;
ALTER TABLE model.build DROP COLUMN date_scheduled;
ALTER TABLE model.build DROP COLUMN date_completed;
ALTER TABLE model.build DROP COLUMN run_by;
ALTER TABLE model.build DROP COLUMN created_by;
ALTER TABLE model.build DROP COLUMN created_at;
ALTER TABLE model.build DROP COLUMN updated_at;

COMMIT;
