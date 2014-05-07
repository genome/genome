-- Deploy timeline_allocation_event_type_values
-- requires: timeline_allocation_event_type

-- It's written this way because these validation values are required for
-- the system to function, and were already in the table, but there was no
-- tracking for how they should get into the table in the first place
BEGIN;

DO $$
    DECLARE vals varchar[] := ARRAY['created', 'purged', 'preserved', 'moved', 'reallocated',
                       'archived', 'unarchived', 'unpreserved', 'finalized',
                        'invalidated', 'strengthened', 'weakened'];
    DECLARE i varchar;
BEGIN
    FOREACH i IN ARRAY vals
    LOOP
        BEGIN
            insert into timeline.allocation_event_type values (i);
        EXCEPTION when unique_violation THEN
            -- do nothing
        END;
    END LOOP;
END;
$$ language 'plpgsql';

COMMIT;
