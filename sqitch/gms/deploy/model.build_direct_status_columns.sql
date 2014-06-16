-- Deploy model.build_direct_status_columns
-- requires: model_build

BEGIN;
    ALTER TABLE model.build ADD COLUMN status TEXT;
    ALTER TABLE model.build ADD COLUMN date_scheduled TIMESTAMP(6) WITHOUT TIME ZONE;
    ALTER TABLE model.build ADD COLUMN date_completed TIMESTAMP(6) WITHOUT TIME ZONE;
    ALTER TABLE model.build ADD COLUMN run_by TEXT;
    ALTER TABLE model.build ADD COLUMN created_by TEXT;
    ALTER TABLE model.build ADD COLUMN created_at TIMESTAMP(6) WITHOUT TIME ZONE;
    ALTER TABLE model.build ADD COLUMN updated_at TIMESTAMP(6) WITHOUT TIME ZONE;

    CREATE INDEX build_status_idx ON model.build USING btree(status);
    CREATE INDEX build_date_scheduled_idx ON model.build USING btree(date_scheduled);
    CREATE INDEX build_date_completed_idx ON model.build USING btree(date_completed);
    CREATE INDEX build_run_by_idx ON model.build USING btree(run_by);
COMMIT;
