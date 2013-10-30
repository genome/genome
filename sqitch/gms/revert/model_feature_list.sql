-- Revert model_feature_list

BEGIN;

DROP TABLE IF EXISTS model.feature_list;

COMMIT;
