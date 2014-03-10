-- Verify model.model.drop_uniq_m_m_idun

BEGIN;

ALTER TABLE model.model ADD CONSTRAINT uniq_m_m_idun UNIQUE (genome_model_id, user_name);

ROLLBACK;
