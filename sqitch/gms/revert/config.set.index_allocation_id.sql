-- Revert config.set.index_allocation_id

BEGIN;

DROP INDEX config.c_s_allocation_id_index;

COMMIT;
