-- Deploy timeline.allocation.absolute_path
-- requires: timeline_allocation

BEGIN;

CREATE INDEX allocation_absolute_path_idx on timeline.allocation using btree (absolute_path);

COMMIT;
