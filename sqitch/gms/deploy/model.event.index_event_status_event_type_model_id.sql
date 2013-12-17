-- Deploy model.event.event_status_event_type_model_id
-- requires: model_event

BEGIN;

CREATE INDEX idx_m_e_es_et_mi on model.event using btree (event_status, event_type, model_id);

COMMIT;
