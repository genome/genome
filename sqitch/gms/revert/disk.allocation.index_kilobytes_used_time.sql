-- Revert disk.allocation.index_kilobytes_used_time

BEGIN;

DROP INDEX disk.d_a_kilobytes_used_time_index;

COMMIT;
