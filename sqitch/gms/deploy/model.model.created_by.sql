-- Deploy model.model.created_by
-- requires: model_model

BEGIN;

ALTER TABLE model.model ADD COLUMN created_by character varying(64);

COMMIT;
