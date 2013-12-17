-- Revert disk.file_summary.index_allocation_id

BEGIN;

DROP INDEX disk.file_summary_allocation_id_idx;

COMMIT;
