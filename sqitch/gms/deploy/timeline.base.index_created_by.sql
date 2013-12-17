-- Deploy timeline.base.created_by
-- requires: timeline_base

BEGIN;

CREATE INDEX base_created_by_idx on timeline.base using btree (created_by);

COMMIT;
