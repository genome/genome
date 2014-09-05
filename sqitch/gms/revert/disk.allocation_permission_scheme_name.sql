-- Revert disk.allocation_permission_scheme_name

BEGIN;

ALTER TABLE disk.allocation_permission_scheme DROP COLUMN name;

COMMIT;
