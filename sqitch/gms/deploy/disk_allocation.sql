-- Deploy disk_allocation
-- requires: disk_schema

BEGIN;

ALTER TABLE disk.allocation DROP COLUMN preserved;
ALTER TABLE disk.allocation ALTER COLUMN archive_after_time SET NOT NULL;
ALTER TABLE disk.allocation ALTER COLUMN status SET NOT NULL;

COMMIT;
