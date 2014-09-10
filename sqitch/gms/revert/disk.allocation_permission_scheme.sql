-- Revert disk.allocation_permission_scheme

BEGIN;

ALTER TABLE disk.allocation DROP COLUMN permission_scheme_id;
DROP TABLE disk.allocation_permission_scheme;

COMMIT;
