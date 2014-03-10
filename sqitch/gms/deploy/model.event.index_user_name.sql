-- Deploy model.event.user_name
-- requires: model_event

BEGIN;

CREATE INDEX event_user_name_index on model.event using btree (user_name);

COMMIT;
