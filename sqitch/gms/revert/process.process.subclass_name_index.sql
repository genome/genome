-- Revert process.process.subclass_name_index

BEGIN;

DROP INDEX process.process_subclass_name_idx;

COMMIT;
