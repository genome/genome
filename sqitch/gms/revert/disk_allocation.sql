-- Deploy disk_allocation
-- requires: disk_schema

BEGIN;

ALTER TABLE disk.allocation ADD COLUMN preserved boolean;
ALTER TABLE disk.allocation ALTER COLUMN archive_after_time DROP NOT NULL;
ALTER TABLE disk.allocation ALTER COLUMN status DROP NOT NULL;

COMMIT;
