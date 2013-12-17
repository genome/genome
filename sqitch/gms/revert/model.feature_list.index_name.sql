-- Revert model.feature_list.index_name

BEGIN;

DROP INDEX model.feature_list_name_index;

COMMIT;
