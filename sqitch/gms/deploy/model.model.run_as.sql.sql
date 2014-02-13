-- Deploy model.model.run_as.sql
-- requires: model_model

BEGIN;

    ALTER TABLE model.model ADD COLUMN run_as VARCHAR(64);

COMMIT;
