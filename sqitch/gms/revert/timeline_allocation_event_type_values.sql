-- Revert timeline_allocation_event_type_values

BEGIN;

DELETE FROM timeline.allocation_event_type where id = 'created';
DELETE FROM timeline.allocation_event_type where id = 'purged';
DELETE FROM timeline.allocation_event_type where id = 'preserved';
DELETE FROM timeline.allocation_event_type where id = 'moved';
DELETE FROM timeline.allocation_event_type where id = 'reallocated';
DELETE FROM timeline.allocation_event_type where id = 'archived';
DELETE FROM timeline.allocation_event_type where id = 'unarchived';
DELETE FROM timeline.allocation_event_type where id = 'unpreserved';
DELETE FROM timeline.allocation_event_type where id = 'finalized';
DELETE FROM timeline.allocation_event_type where id = 'invalidated';
DELETE FROM timeline.allocation_event_type where id = 'strengthened';
DELETE FROM timeline.allocation_event_type where id = 'weakened';

COMMIT;
