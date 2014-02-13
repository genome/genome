-- Revert model.model.created_by

BEGIN;

ALTER TABLE model.model DROP COLUMN IF EXISTS created_by;

COMMIT;
