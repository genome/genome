-- Revert timeline_allocation_event_type_values

BEGIN;

select 1/count(*) from timeline.allocation_event_type where id = 'created';
select 1/count(*) from timeline.allocation_event_type where id = 'purged';
select 1/count(*) from timeline.allocation_event_type where id = 'preserved';
select 1/count(*) from timeline.allocation_event_type where id = 'moved';
select 1/count(*) from timeline.allocation_event_type where id = 'reallocated';
select 1/count(*) from timeline.allocation_event_type where id = 'archived';
select 1/count(*) from timeline.allocation_event_type where id = 'unarchived';
select 1/count(*) from timeline.allocation_event_type where id = 'unpreserved';
select 1/count(*) from timeline.allocation_event_type where id = 'finalized';
select 1/count(*) from timeline.allocation_event_type where id = 'invalidated';
select 1/count(*) from timeline.allocation_event_type where id = 'strengthened';
select 1/count(*) from timeline.allocation_event_type where id = 'weakened';

COMMIT;
