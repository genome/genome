-- Revert model.event.index_prior_event_id

BEGIN;

DROP INDEX model.event_prior_event_index;

COMMIT;
