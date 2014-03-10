-- Deploy model.event.build_id_date_completed_event_type
-- requires: model_event

BEGIN;

CREATE INDEX idx_m_e_bi_dc_et on model.event using btree (build_id, date_completed, event_type);

COMMIT;
