-- Revert result.software_result.index_lookup_hash_class_name

BEGIN;

DROP INDEX result.lookup_hash_class_name_idx;

COMMIT;
