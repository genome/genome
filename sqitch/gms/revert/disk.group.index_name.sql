-- Revert disk.group.index_name

BEGIN;

DROP INDEX disk.group_name_index;

COMMIT;
