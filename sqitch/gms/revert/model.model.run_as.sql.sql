-- Revert model.model.run_as.sql

BEGIN;

    ALTER TABLE model.model DROP COLUMN IF EXISTS run_as;

COMMIT;
