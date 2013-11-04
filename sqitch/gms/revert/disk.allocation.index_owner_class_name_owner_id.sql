-- Revert disk.allocation.index_owner_class_name_owner_id

BEGIN;

DROP INDEX disk.allocation_owner_class_id_index;

COMMIT;
