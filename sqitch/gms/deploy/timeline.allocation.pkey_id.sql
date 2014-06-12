-- Deploy timeline.allocation.pkey_id
-- requires: timeline_allocation

CREATE UNIQUE INDEX CONCURRENTLY allocation_pkey_idx ON timeline.allocation (id);
BEGIN;
    ALTER TABLE timeline.allocation ADD CONSTRAINT allocation_pkey PRIMARY KEY USING INDEX allocation_pkey_idx;
COMMIT;
