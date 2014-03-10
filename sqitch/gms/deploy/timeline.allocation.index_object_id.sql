-- Deploy timeline.allocation.object_id
-- requires: timeline_allocation

BEGIN;

CREATE INDEX allocation_object_id_idx on timeline.allocation using btree (object_id);

COMMIT;
