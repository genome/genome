-- Deploy model.model.drop_uniq_m_m_idun
-- requires: model_model

BEGIN;

ALTER TABLE model.model DROP CONSTRAINT uniq_m_m_idun;

COMMIT;
