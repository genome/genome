-- Deploy timeline.allocation.name
-- requires: timeline_allocation

BEGIN;

CREATE INDEX allocation_name_idx on timeline.allocation using btree (name);

COMMIT;
