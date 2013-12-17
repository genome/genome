-- Deploy workflow_plan
-- requires: workflow_schema

BEGIN;

CREATE TABLE workflow.plan (
    workflow_plan_id character varying(255) NOT NULL,
    xml xml,
    CONSTRAINT plan_pk PRIMARY KEY (workflow_plan_id)
)
WITH (autovacuum_enabled=false, toast.autovacuum_enabled=false);

COMMIT;
