-- Verify model.event.index_prior_event_id

BEGIN;

SELECT 1/count(*) FROM pg_class WHERE relkind = 'i' and relname = 'event_prior_event_index';

ROLLBACK;
