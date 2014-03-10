-- Deploy model.model.build_requested
-- requires: model_model

BEGIN;

CREATE INDEX m_m_build_requested_index on model.model using btree (build_requested);

COMMIT;
