-- Deploy model.event.prior_event_id
-- requires: model_event

BEGIN;

CREATE INDEX event_prior_event_index on model.event using btree (prior_event_id);

COMMIT;
