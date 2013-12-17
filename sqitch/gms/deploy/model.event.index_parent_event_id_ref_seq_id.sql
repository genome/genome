-- Deploy model.event.parent_event_id_ref_seq_id
-- requires: model_event

BEGIN;

CREATE INDEX idx_m_e_pei_rsi on model.event using btree (parent_event_id, ref_seq_id);

COMMIT;
