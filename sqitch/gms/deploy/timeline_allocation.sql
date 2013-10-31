-- Deploy timeline_allocation
-- requires: timeline_base
-- requires: disk_allocation
-- requires: timeline_allocation_event_type

BEGIN;

CREATE TABLE IF NOT EXISTS timeline.allocation (
    kilobytes_requested character varying NOT NULL,
    absolute_path character varying NOT NULL,
    status character varying NOT NULL,
    CONSTRAINT allocation_event_allocation_fk FOREIGN KEY (object_id) REFERENCES disk.allocation(id) MATCH FULL,
    CONSTRAINT allocation_event_typ_fk FOREIGN KEY (name) REFERENCES timeline.allocation_event_type(id) MATCH FULL
)
INHERITS (timeline.base);

COMMIT;
