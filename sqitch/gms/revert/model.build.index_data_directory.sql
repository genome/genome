-- Revert model.build.index_data_directory

BEGIN;

DROP INDEX model.build_directory_index;

COMMIT;
