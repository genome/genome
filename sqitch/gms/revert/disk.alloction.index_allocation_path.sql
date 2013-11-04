-- Revert disk.alloction.index_allocation_path

BEGIN;

DROP INDEX disk.idx_allocation_path_like;

COMMIT;
